#include <array>
#include <chrono>
// #include <cstddef>
#include <string>
// #include <type_traits>
#include <vector>

#include "analytic/0q3gA-analytic.h"
#include "chsums/0q3gA.h"
#include "ngluon2/Model.h"
#include "ngluon2/Mom.h"
// #include "ngluon2/refine.h"
#include "tools/PhaseSpace.h"

template <typename T>
class Test {
private:
    static constexpr int mul { 5 };
    const double MuR2 { std::pow(91.118, 2) };

    Amp0q3gAA<T> num;
    Amp0q3gAA_a<T> ana;

    const std::vector<MOM<T>> genMom(int rseed) const
    {
        std::array<Flavour<double>, mul> massless;
        massless.fill(StandardModel::G());

        // args: multiplicity, flavour array, rseed, sqrts
        PhaseSpace<double> psp(mul, massless.data(), rseed, 1.);
        std::vector<MOM<double>> mom_double { psp.getPSpoint() };

#ifdef DEBUG
        psp.showPSpoint();
#endif // DEBUG

        // std::vector<MOM<T>> mom(mul);
        // refine(mom_double, mom);
        // return mom;

        return mom_double;
    }

public:
    const std::string genHelStr(const std::array<int, mul>& helInt) const
    {
        std::string helStr;
        for (int hel : helInt) {
            helStr.append(hel == 1 ? "+" : "-");
        }
        return helStr;
    }

    constexpr int get_mul() const
    {
        return this->mul;
    }

    std::chrono::nanoseconds compare(std::array<int, mul> hel, const bool v = true)
    {
        if (v) {
            std::cout << this->genHelStr(hel) << '\n';
        }

        std::chrono::high_resolution_clock::time_point t0s { std::chrono::high_resolution_clock::now() };
        const T numsol { this->num.virtsq(hel.data()).get0().real() };
        std::chrono::high_resolution_clock::time_point t0e { std::chrono::high_resolution_clock::now() };
        std::chrono::nanoseconds d0 { std::chrono::duration_cast<std::chrono::nanoseconds>(t0e - t0s).count() };
        if (v) {
            std::cout << "|amp_n|^2 " << numsol << '\n';
            std::cout << "t_n       " << std::chrono::duration_cast<std::chrono::milliseconds>(d0).count() << "ms\n";
        }

        std::chrono::high_resolution_clock::time_point t1s { std::chrono::high_resolution_clock::now() };
        const T anasol { this->ana.virtsq(hel.data()).get0().real() };
        std::chrono::high_resolution_clock::time_point t1e { std::chrono::high_resolution_clock::now() };
        std::chrono::nanoseconds d1 { std::chrono::duration_cast<std::chrono::nanoseconds>(t1e - t1s).count() };
        if (v) {
            std::cout << "|amp_a|^2 " << anasol << '\n';
            std::cout << "t_a       " << std::chrono::duration_cast<std::chrono::microseconds>(d1).count() << "\u03BCs\n";
        }

        const long int speedup { d0 / d1 };
        if (v) {
            std::cout << "ratio     " << anasol / numsol << '\n';
            std::cout << "speedup   " << speedup << '\n';
        }

        return d1;
    }

    void primary_helicities()
    {

        // // +++++
        // std::array<int, this->mul> hel3 { 1, 1, 1, 1, 1 };
        // this->compare(hel3);

        // // -++++
        // std::array<int, this->mul> hel4 { -1, 1, 1, 1, 1 };
        // this->compare(hel4);

        // --+++
       //std::array<int, this->mul> hel1 { -1, -1, 1, 1, 1 };
       // this->compare(hel1);

        // +++--
        std::array<int, this->mul> hel2 { 1, 1, 1, -1, -1 };
        this->compare(hel2);

        // // +++-+
        // std::array<int, this->mul> hel5 { 1, 1, 1, -1, 1 };
        // this->compare(hel5);

        // +--+-
       // std::array<int, this->mul> hel6 { 1, -1, -1, 1, -1 };
       // this->compare(hel6);
    };

    void all_helicities()
    {
        for (int i0 { -1 }; i0 < 2; i0 += 2) {
            for (int i1 { -1 }; i1 < 2; i1 += 2) {
                for (int i2 { -1 }; i2 < 2; i2 += 2) {
                    for (int i3 { -1 }; i3 < 2; i3 += 2) {
                        for (int i4 { -1 }; i4 < 2; i4 += 2) {
                            const std::array<int, mul> hel { { i4, i3, i2, i1, i0 } };
                            this->compare(hel);
                        }
                    }
                }
            }
        }
    }

    Test(int rseed)
        : num(StandardModel::Ax(StandardModel::IL(), StandardModel::IL().C()), 1.)
        , ana(StandardModel::Ax(StandardModel::IL(), StandardModel::IL().C()), 1.)
    {
        const std::vector<MOM<T>> mom { genMom(rseed) };

        this->num.setMomenta(mom);
        this->num.setMuR2(this->MuR2);

        this->ana.setMomenta(mom);
        this->ana.setMuR2(this->MuR2);
    }
};

int main(int argc, char** argv)
{
    std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
    std::cout.precision(16);

    if (argc != 1) {
        std::cout << "Warning: " << argv[0] << " doesn't accept command line arguments.\n";
    }

    const int numPSPs { 1 };
    const int rseed1 { 1 };
    // std::chrono::nanoseconds total { 0 };
    for (int i { rseed1 }; i < numPSPs + rseed1; ++i) {
        Test<double> test(i);
        test.primary_helicities();
        // test.all_helicities();

        // std::array<int, test.get_mul()> hel { -1, -1, 1, 1, 1 };
        // total += test.compare(hel, true);
    }
    // std::cout << "<t_a>     " << std::chrono::duration_cast<std::chrono::microseconds>(total / numPSPs).count() << "\u03BCs\n";

    return 0;
}

