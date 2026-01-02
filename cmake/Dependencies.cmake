include(FetchContent)

FetchContent_Declare(
    core-math
    GIT_REPOSITORY ssh://git@github.com/marquesasarbitrage/core-math.git
    GIT_TAG main  # or a release tag
)
FetchContent_MakeAvailable(core-math)

FetchContent_Declare(
    core-datetime
    GIT_REPOSITORY ssh://git@github.com/marquesasarbitrage/core-datetime.git
    GIT_TAG main  # or a release tag
)
FetchContent_MakeAvailable(core-datetime)