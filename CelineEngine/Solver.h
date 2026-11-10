#pragma once

#include <vector>

namespace CelineEngine
{
    //==========================================================================
    /**
        A small dense LU solver with partial pivoting, sized once and then
        reused forever without allocating.

        Circuit solves G·V = I every sample, and once per Newton iteration when
        nonlinear devices are involved. All storage is allocated in resize(),
        from Circuit::prepare() on the message thread; every method afterwards is
        allocation-free.

        Two phases, so the linear case skips most of the work:

            std::copy (stampedMatrix, ..., solver.data());
            solver.factorise();      // O(n^3), only when values change
            solver.solveInPlace (b); // O(n^2), every sample

        Double precision on purpose: a circuit mixes conductances over many
        orders of magnitude -- gmin at 1e-12 S, a 1M resistor at 1e-6 S, a
        conducting diode at ~1e-1 S -- and float's seven digits are not enough
        headroom once diodes are in the matrix. Only the solve is double; the
        audio in and out stays float.
    */
    class DenseSolver
    {
       public:
        /** Allocates for an n x n system. Call from prepare(), never from the audio thread. */
        void resize(int n);

        int size() const noexcept { return dimension; }

        /** Raw matrix storage. Callers stamp a full system elsewhere and bulk-copy
            it in, which is why there's no element accessor here: every caller
            already has its own matrix to build, and a second stamping path would
            be one more thing to keep in step. */
        double* data() noexcept { return matrix.data(); }
        const double* data() const noexcept { return matrix.data(); }

        /** In-place LU decomposition with partial pivoting.
            Returns false if the matrix is singular, in which case the stored
            factorisation must not be used. */
        bool factorise() noexcept;

        /** Solves for the right-hand side using the factorisation from the last
            successful factorise() call. `rhs` (n elements) is overwritten with
            the solution. */
        void solveInPlace(double* rhs) const noexcept;

        /** Factorises and solves in one pass, for a matrix that is used for a
            single solve and then rebuilt -- a Newton iteration's Jacobian.

            The right-hand side rides along with the elimination, so the forward
            substitution happens inside the row operations rather than after
            them, and none of what solveInPlace() needs for later solves is
            built. Which is also why solveInPlace() must not be called after
            this: the factorisation it would use is not there.

            `rhs` (n elements) is overwritten with the solution. Returns false
            if the matrix is singular. */
        bool factoriseAndSolve(double* rhs) noexcept;

        /** Solves the n x n row-major system a x = rhs in place by Gaussian
            elimination with partial pivoting, rhs becoming x and a destroyed.
            1/U(k,k) goes into `invDiagonal` when given, otherwise onto a's own
            diagonal, which needs no storage of its own -- for a system whose size
            changes from call to call, where a sized DenseSolver would allocate.

            The one statement of the elimination: factoriseAndSolve() and the
            reduced Newton step both call it, and a test holds them to the bit.
            False if a pivot vanishes. */
        static bool eliminateAndSolve(double* a, double* rhs, int n,
                                      double* invDiagonal = nullptr) noexcept;

       private:
        /** The row at or below k with the largest entry in column k. */
        static int pivotRowFor(const double* a, int n, int k) noexcept;

        /** A pivot smaller than this means the matrix is singular. */
        static constexpr double singularPivot = 1.0e-300;

        int dimension = 0;
        std::vector<double> matrix;      // n x n, row-major; holds L and U after factorise()
        std::vector<int> rowSwaps;       // pivot row chosen at each elimination step
        std::vector<double> invDiagonal; // 1/U(k,k), so back substitution multiplies
        std::vector<double> transposed;  // the factors again, column by column, for solveInPlace()
    };
} // namespace CelineEngine
