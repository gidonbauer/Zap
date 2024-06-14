#include <algorithm>
#include <chrono>

#include <csignal>
#include <sys/wait.h>
#include <unistd.h>

#include "Igor.hpp"
#include "ReadMatrixInc.hpp"
#include "Renderer/Canvas.hpp"

namespace Rd = Zap::Renderer;

constexpr auto READ_END  = 0UZ;
constexpr auto WRITE_END = 1UZ;

auto start_ffmpeg(size_t width,
                  size_t height,
                  std::string_view output_file) noexcept -> std::pair<pid_t, pid_t> {
  std::array<int, 2> pipefd{};
  if (pipe(pipefd.data()) == -1) {
    Igor::Panic("Opening pipe failed: {}", std::strerror(errno));
  }

  const pid_t child = fork();
  if (child == -1) {
    Igor::Panic("Could not fork, no child process created: {}", std::strerror(errno));
  }

  if (child == 0) {
    if (dup2(pipefd[READ_END], STDIN_FILENO) == -1) {
      Igor::Panic("Could not attach read end of pipe in child process: {}", std::strerror(errno));
    }
    if (close(pipefd[WRITE_END]) == -1) {
      Igor::Panic("Could not close write end of pipe in child process: {}", std::strerror(errno));
    }

    Igor::Info("Start FFmpeg...");

    const auto resolution = std::format("{}x{}", width, height);
    const auto blocksize  = std::to_string(width * height * sizeof(Rd::Canvas::PixelType));
    static_assert(sizeof(Rd::Canvas::PixelType) == 24UZ / 8UZ,
                  "FFmpeg assumes that the pixels are in the rgb24 format, but the size of the "
                  "pixel datatype of the canvas is too large.");
    // clang-format off
    int ret = execlp("ffmpeg", // NOLINT
        "ffmpeg",
        // "-loglevel", "verbose",
        "-y",
        "-f", "rawvideo",
        "-pix_fmt", "rgb24",
        "-s", resolution.c_str(),
        "-r", "60",
        "-blocksize", blocksize.c_str(),
        "-i", "-",
        "-c:v", "libx264",
        "-vb", "2500k",
        "-pix_fmt", "yuv420p",
        output_file.data(),
        static_cast<const char*>(nullptr));
    // clang-format on

    if (ret == -1) {
      Igor::Panic("An error occured in `execlp`: {}", std::strerror(errno));
    }
  } else {
    if (close(pipefd[READ_END]) == -1) {
      Igor::Panic("Could not close read end of pipe in parent process: {}", std::strerror(errno));
    }

    return {child, pipefd[WRITE_END]};
  }
  Igor::Panic("Unreachable");
  std::unreachable();
}

auto end_ffmpeg(pid_t child, pid_t stream) noexcept -> bool {
  if (close(stream) == -1) {
    Igor::Warn("Could not close write end of pipe in parent process: {}", std::strerror(errno));
    return false;
  }

  int status = 0;
  if (waitpid(child, &status, 0) == -1) {
    Igor::Warn("Could not wait for child process to finish: {}", std::strerror(errno));
    return false;
  }
  Igor::Info("FFmpeg process finished with status {}.", status);

  return true;
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  using Float = double;

  if (argc < 4) {
    Igor::Warn("Usage: {} <u input file> <t input file> <output_file>", *argv);
    return 1;
  }
  const auto u_input_file = argv[1];  // NOLINT
  const auto t_input_file = argv[2];  // NOLINT
  const auto output_file  = argv[3];  // NOLINT

  Zap::IncMatrixReader<Float> u_reader{u_input_file};
  Zap::IncMatrixReader<Float> t_reader{t_input_file};
  if (u_reader.num_elem() != t_reader.num_elem()) {
    Igor::Warn("Incompatible number of elements in u ({}) and t ({}).", u_input_file, t_input_file);
    return 1;
  }

  constexpr size_t TEXT_HEIGHT      = 100UZ;
  constexpr size_t MIN_CANVAS_WIDTH = 600UZ;
  Rd::Canvas canvas(std::max(MIN_CANVAS_WIDTH, static_cast<size_t>(u_reader.cols())),
                    static_cast<size_t>(u_reader.rows()) + TEXT_HEIGHT);
  Rd::Box text_box{
      .col    = 0,
      .row    = 0,
      .width  = canvas.width(),
      .height = TEXT_HEIGHT,
  };

  const size_t avail_padding = canvas.width() - static_cast<size_t>(u_reader.cols());
  Rd::Box graph_box{
      .col    = avail_padding / 2UZ,
      .row    = text_box.height,
      .width  = static_cast<size_t>(u_reader.cols()),
      .height = static_cast<size_t>(u_reader.rows()),
  };
  constexpr std::string_view font_file = "../assets/LinLibertine_R.ttf";
  if (!canvas.load_font(font_file)) {
    Igor::Warn("Could not load font from file `{}`.", font_file);
    return 1;
  }

  const auto [ffmpeg_pid, ffmpeg_stream] =
      start_ffmpeg(canvas.width(), canvas.height(), output_file);

  const auto t_begin = std::chrono::high_resolution_clock::now();
  while (u_reader.read_next<false>() && t_reader.read_next<false>()) {
    canvas.clear({});

    if (const std::string str = std::format("t = {:.6f}", t_reader(0, 0));
        !canvas.draw_text(str, text_box, true)) {
      Igor::Warn("Could not render string `{}` to canvas.", str);
      goto error_after_ffmpeg;  // NOLINT
    }

    const auto& data      = u_reader.data();
    const auto [min, max] = std::minmax_element(std::cbegin(data), std::cend(data));
    if (!canvas.draw_buffer(
            data,
            graph_box,
            [&](Float v) { return Rd::float_to_rgb<Rd::Float2RGB::COLORMAP>(v, *min, *max); },
            u_reader.is_row_major(),
            true)) {
      Igor::Warn("Could not draw graph to canvas.");
      goto error_after_ffmpeg;  // NOLINT
    }

    if (!canvas.to_raw_stream(ffmpeg_stream)) {
      Igor::Warn("Could not write canvas to FFmpeg.");
      goto error_after_ffmpeg;  // NOLINT
    }
  }

  if (!end_ffmpeg(ffmpeg_pid, ffmpeg_stream)) {
    Igor::Warn("Could not properly end FFmpeg.");
    goto error_after_ffmpeg;  // NOLINT
  }
  {
    const auto t_dur =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_begin);
    Igor::Info("Rendering took {}.", t_dur);
  }

  return 0;
error_after_ffmpeg:
  if (kill(ffmpeg_pid, SIGKILL) == -1) {
    Igor::Warn("Could not kill FFmpeg: {}", std::strerror(errno));
  }
  return 1;
}
