#ifndef ZAP_RENDERER_FFMPEG_HPP_
#define ZAP_RENDERER_FFMPEG_HPP_

#include <array>
#include <cerrno>
#include <cstring>
#include <string_view>

#include <csignal>
#include <sys/wait.h>
#include <unistd.h>

#include "Igor.hpp"
#include "Renderer/Canvas.hpp"

namespace Zap::Renderer {

class FFmpeg {
  enum : size_t { READ_END = 0UZ, WRITE_END = 1UZ };
  pid_t m_child{};
  pid_t m_stream{};

 public:
  [[nodiscard]] FFmpeg(size_t width,
                       size_t height,
                       const std::string& output_file,
                       size_t fps = 60) noexcept {
    std::array<int, 2> pipefd{};
    if (pipe(pipefd.data()) == -1) { Igor::Panic("Opening pipe failed: {}", std::strerror(errno)); }

    m_child = fork();
    if (m_child == -1) {
      Igor::Panic("Could not fork, no child process created: {}", std::strerror(errno));
    }

    if (m_child == 0) {
      if (dup2(pipefd[READ_END], STDIN_FILENO) == -1) {
        Igor::Panic("Could not attach read end of pipe in child process: {}", std::strerror(errno));
      }
      if (close(pipefd[WRITE_END]) == -1) {
        Igor::Panic("Could not close write end of pipe in child process: {}", std::strerror(errno));
      }

      Igor::Info("Start FFmpeg...");

      const auto resolution = std::format("{}x{}", width, height);
      const auto blocksize  = std::to_string(width * height * sizeof(Canvas::PixelType));
      static_assert(sizeof(Canvas::PixelType) == 24UZ / 8UZ,
                    "FFmpeg assumes that the pixels are in the rgb24 format, but the size of the "
                    "pixel datatype of the canvas is too large.");

      assert(fps > 0);
      const auto FPS = std::to_string(fps);

      // clang-format off
    int ret = execlp("ffmpeg", // NOLINT
        "ffmpeg",
        // "-loglevel", "verbose",
        "-y",
        "-f", "rawvideo",
        "-pix_fmt", "rgb24",
        "-s", resolution.c_str(),
        "-r", FPS.c_str(),
        "-blocksize", blocksize.c_str(),
        "-probesize", "10M",
        "-i", "-",
        "-c:v", "libx264",
        "-b:v", "3500k",
        "-pix_fmt", "yuv420p",
        output_file.c_str(),
        static_cast<const char*>(nullptr));
      // clang-format on

      if (ret == -1) { Igor::Panic("An error occured in `execlp`: {}", std::strerror(errno)); }
    } else {
      if (close(pipefd[READ_END]) == -1) {
        Igor::Panic("Could not close read end of pipe in parent process: {}", std::strerror(errno));
      }

      m_stream = pipefd[WRITE_END];
    }
  }

  FFmpeg(const FFmpeg& other) noexcept                    = delete;
  FFmpeg(FFmpeg&& other) noexcept                         = delete;
  auto operator=(const FFmpeg& other) noexcept -> FFmpeg& = delete;
  auto operator=(FFmpeg&& other) noexcept -> FFmpeg&      = delete;

  ~FFmpeg() noexcept {
    bool encountered_error = false;
    if (close(m_stream) == -1) {
      Igor::Warn("Could not close write end of pipe in parent process: {}", std::strerror(errno));
      encountered_error = true;
    }

    int status = 0;
    if (waitpid(m_child, &status, 0) == -1) {
      Igor::Warn("Could not wait for child process to finish: {}", std::strerror(errno));
      encountered_error = true;
    }

    if (encountered_error) {
      if (kill(m_child, SIGKILL) == -1) {
        Igor::Warn("Could not kill FFmpeg: {}", std::strerror(errno));
      }
      Igor::Info("FFmpeg process has been terminated.");
    } else {
      Igor::Info("FFmpeg process finished with status {}.", status);
    }
  }

  [[nodiscard]] constexpr auto stream() const noexcept -> pid_t { return m_stream; }
};

}  // namespace Zap::Renderer

#endif  // ZAP_RENDERER_FFMPEG_HPP_
