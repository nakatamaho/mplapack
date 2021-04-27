#include <fem.hpp> // Fortran EMulation library of fable module

using namespace fem::major_types;
using fem::common;

int main(int argc, char const* argv[])
{
  common cmn(argc, argv);
  common_write write(cmn);
  write(6, star), "HELLO, WORLD!";
}

