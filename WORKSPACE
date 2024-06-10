workspace(name = "highs")

load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository", "new_git_repository")

# Needed by Protobuf
git_repository(
    name = "rules_python",
    tag = "0.31.0",
    remote = "https://github.com/bazelbuild/rules_python.git",
)

load("@rules_python//python:repositories.bzl", "py_repositories")
py_repositories()

# Protobuf
git_repository(
    name = "com_google_protobuf",
    tag = "v26.1",
    remote = "https://github.com/protocolbuffers/protobuf.git",
)

# Load common dependencies.
load("@com_google_protobuf//:protobuf_deps.bzl", "protobuf_deps")
protobuf_deps()

# ZLIB
new_git_repository(
    name = "zlib",
    build_file = "@com_google_protobuf//:third_party/zlib.BUILD",
    tag = "v1.2.11",
    remote = "https://github.com/madler/zlib.git",
)

