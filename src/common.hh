#ifndef __common_hh__
#define __common_hh__

#include <vector>
#include <string>

void check_boundary_vector(std::string name, int n, const std::vector<double>& v);

void convolve_same_size(int size, const double* src0, const double* src1, double* dest);

template<class T>
class DoubleBuffer {
    public:
        DoubleBuffer(int n, T value);
        std::vector<T>& get_src()  { return buf0_is_src ? buf0 : buf1; };
        std::vector<T>& get_dest() { return buf0_is_src ? buf1 : buf0; };
        void flip() { buf0_is_src = !buf0_is_src; };

    private:
        std::vector<T> buf0;
        std::vector<T> buf1;
        bool buf0_is_src;
};

template<class T>
DoubleBuffer<T>::DoubleBuffer(int n, T value) :
    buf0(n, value), buf1(n, value), buf0_is_src(true)
{
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << std::string("[");
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << std::string("\b\b]");
  }
  return out;
}

#endif

