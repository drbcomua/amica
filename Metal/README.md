## Install dependencies (MacOS)

```shell
git clone https://github.com/metal-cpp/metal-cpp.git
```

## Compile

```shell
clang++ amicable_metal.cpp -o amicable_finder_metal \
        -std=c++17 -O3 -pthread \
        -I./metal-cpp \
        -framework Foundation -framework Metal -framework QuartzCore
```

## Run

```shell
# Example: Search for pairs up to 10,000,000
./amicable_finder_metal 10000000
```

## Expected result

```shell
```