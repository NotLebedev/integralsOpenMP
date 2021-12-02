#include <iostream>
#include <mpi.h>

#include "types.h"
#include "functions/arcsin.h"
#include "functions/heaviside_step.h"
#include "functions/exp.h"
#include "messaging.h"

/**
 * Compound function with heavy calculations
 */
data_t func_big(data_t x) {
    return exp(x) + heaviside_step(x) + arcsin(x);
}

bool invalidate = false;

class Result {
public:
    void timer_start() {
        time = MPI_Wtime();
    }

    void timer_end() {
        time = MPI_Wtime(); - time;
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

data_t run(const Partition x_, const func_type func) {
    data_t result = 0.0;

    size_t n = x_.get_n();
    for (size_t i = 1; i < n - 1; i++) {
        result += func(x_(i));
    }

    result += (func(x_(n)) + func(x_(0))) / 2;
    result *= x_.get_delta();

    return result;
}

Result run_full(const size_t n, const data_t a, const data_t b, const func_type func) {
    Result res{};
    res.timer_start();

    Partition full{a, b, n};

    int num_process;
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);

    size_t partition_step = n / num_process;
    for (size_t i = 0, j = 0; j < num_process; j++) {
        if (j < n % num_process) {
            Partition p{full(i), full(i + partition_step + 1), partition_step + 1};
            //printf("Partition %lu %lu %lu\n", i, i + partition_step + 1, partition_step + 1);
            i += partition_step + 1;
            master_send_job(p, (int) j);
        } else if (partition_step != 0) {
            Partition p{full(i), full(i + partition_step), partition_step};
            //printf("Partition %lu %lu %lu\n", i, i + partition_step, partition_step);
            i += partition_step;
            master_send_job(p, (int) j);
        }
    }

    Partition p{full(0), full(partition_step + (n % num_process > 0)), partition_step + (n % num_process > 0)};
    data_t r = run(p, func);
    r = reduce(r);

    res.set_result(r);
    res.timer_end();

    return res;
}

void benchmark(const size_t n, const data_t a, const data_t b, const func_type func) {

#ifndef WARMUP_ROUNDS_CNT
#define WARMUP_ROUNDS_CNT 3
#endif
    // Warmup round
    for (size_t i = 0; i < WARMUP_ROUNDS_CNT; i++)
        run_full(n, a, b, func);

#ifndef ROUNDS_CNT
#define ROUNDS_CNT 5
#endif
    // Actual round
    double avg = 0.0;
    for (int i = 0; i < ROUNDS_CNT; i++) {
        // Invalidate coefficient table
        invalidate = true;
        arcsin(0);
        exp(0);
        heaviside_step(0);
        invalidate = false;

        Result res = run_full(n, a, b, func);
        avg  = (avg * i + res.get_runtime()) / (i + 1);
    }

    std::cout << "Average runtime: " << avg << " seconds" << std::endl;
}

void actual(const size_t n, const data_t a, const data_t b, const func_type func) {
    Result res = run_full(n, a, b, func);

    std::cout << "Runtime: " << res.get_runtime() << " seconds. ";
    std::cout << "Result: " << res.get_result() << std::endl;
}

void run_slave() {
    printf("Slave entered\n");
    fflush(stdout);

    while (true) {
        try {
            const Partition p = slave_receive_job();
            const func_type func = func_big;
            data_t res = run(p, func);
            reduce(res);
        } catch (std::exception &e) {
            return;
        }
    }
}

void run_master(int argc, char *argv[]) {

    printf("Master entered\n");
    fflush(stdout);

    size_t n = 1000000;
    if (argc == 2) {
        char *end;
        size_t cnt = strtoull(argv[1], &end, 10);
        if (end != argv[1])
            n = cnt;
    }
    const data_t a = -1.0;
    const data_t b = 1.0;
    const func_type func = func_big;

#ifdef BENCHMARK
    benchmark(n, a, b, func);
#else
    actual(n, a, b, func);
#endif
    int num_process;
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    for (int i = 1; i < num_process; ++i) {
        master_send_terminate(i);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc,&argv);
    int mpi_proc_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_id);
    if (mpi_proc_id != 0)
        run_slave();
    else
        run_master(argc, argv);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}