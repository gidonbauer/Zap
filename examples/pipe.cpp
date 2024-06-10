#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <thread>
using namespace std::chrono_literals;

#include <csignal>
#include <sys/wait.h>
#include <unistd.h>

#include "Igor.hpp"

constexpr auto READ_END  = 0UZ;
constexpr auto WRITE_END = 1UZ;

auto main() -> int {
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
    close(pipefd[WRITE_END]);

    // NOLINTNEXTLINE
    int ret = execlp("/bin/cat", "cat", static_cast<const char*>(nullptr));

    if (ret == -1) {
      Igor::Panic("An error occured in `execlp`: {}", std::strerror(errno));
    }
  } else {
    close(pipefd[READ_END]);

    constexpr std::string_view msg = "Hello world.\n";
    for (int i = 0; i < 5; ++i) {
      write(pipefd[WRITE_END], msg.data(), msg.size());
    }

    close(pipefd[WRITE_END]);

    int status = 0;
    if (waitpid(child, &status, 0) == -1) {
      Igor::Panic("Could not wait for child process to finish: {}", std::strerror(errno));
    }

    Igor::Info("Child process exited with status {}", status);
    Igor::Info("Parent exited successfully.");
  }
}
