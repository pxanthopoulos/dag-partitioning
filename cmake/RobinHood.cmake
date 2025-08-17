# Robin Hood Hashing - Fast & memory efficient hashtable based on robin hood hashing

set(ROBIN_HOOD_VERSION "3.11.5")

FetchContent_Declare(
    robin_hood
    GIT_REPOSITORY https://github.com/martinus/robin-hood-hashing.git
    GIT_TAG ${ROBIN_HOOD_VERSION}
)

FetchContent_MakeAvailable(robin_hood)