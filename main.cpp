#include <cmath>
#include <iostream>
#include <omp.h>

#include "types.h"
#include "functions/arcsin.h"

bool invalidate = false;

class Partition {
public:
    Partition(data_t a, data_t b, size_t n_steps) : a_{a}, b_{b}, n_steps_{n_steps},
                                                    delta_{(b - a) / (data_t)n_steps} {}
    data_t operator() (size_t i) const {
        return a_ + ((double) i) * delta_;
    }

    data_t get_delta() const {
        return delta_;
    }

private:
    data_t a_;
    data_t b_;
    size_t n_steps_;
    data_t delta_;
};

class Result {
public:
    void timer_start() {
        time = omp_get_wtime();
    }

    void timer_end() {
        time = omp_get_wtime() - time;
    }

    void set_result(data_t res) noexcept {
        result = res;
    }

    data_t get_result() const noexcept {
        return result;
    }

    double get_runtime() const noexcept {
        return time;
    }
private:
    double time;
    data_t result;
};

Result run(const size_t n, const data_t a, const data_t b, const func_type func) {
    Result res{};
    res.timer_start();
    data_t result = 0.0;

    const Partition x_(a, b, n);
#pragma omp parallel for shared(func, n, x_), reduction(+:result), default(none)
    for (size_t i = 1; i < n - 1; i++) {
        result += func(x_(i));
    }

    result += (func(x_(n)) + func(x_(0))) / 2;
    result *= x_.get_delta();

    res.timer_end();
    res.set_result(result);

    return res;
}

void benchmark(const size_t n, const data_t a, const data_t b, const func_type func) {

#ifndef WARMUP_ROUNDS_CNT
#define WARMUP_ROUNDS_CNT 10
#endif
    // Warmup round
    for (size_t i = 0; i < WARMUP_ROUNDS_CNT; i++)
        run(n, a, b, func);

#ifndef ROUNDS_CNT
#define ROUNDS_CNT 10
#endif
    // Actual round
    double avg = 0.0;
    for (int i = 0; i < ROUNDS_CNT; i++) {
        // Invalidate coefficient table
        invalidate = true;
        arcsin(0);
        invalidate = false;

        Result res = run(n, a, b, func);
        avg  = (avg * i + res.get_runtime()) / (i + 1);
    }

    std::cout << "Average runtime: " << avg << " seconds" << std::endl;
}

void actual(const size_t n, const data_t a, const data_t b, const func_type func) {
    Result res = run(n, a, b, func);

    std::cout << "Runtime: " << res.get_runtime() << " seconds. ";
    std::cout << "Result: " << res.get_result() << std::endl;
}

int main() {
    const size_t n = 100000;
    const data_t a = 0;
    const data_t b = 1.0;
    const func_type func = arcsin;

#ifdef BENCHMARK
    benchmark(n, a, b, func);
#else
    actual(n, a, b, func);
#endif

}