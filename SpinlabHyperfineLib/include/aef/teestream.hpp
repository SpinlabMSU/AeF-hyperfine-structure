/*
    aef/teestream.hpp -- implements a stream that tees its input to two outputs.

    This file comes from https://stackoverflow.com/a/27216591 and hence is
    under the CC BY-SA 3.0. See https://creativecommons.org/licenses/by-sa/3.0/
    for details.
*/
#pragma once
// from https://stackoverflow.com/a/27216591

#include <streambuf>
#include <ostream>
#include <filesystem>
namespace teestream {
    template <class charT, class traits = std::char_traits<charT> >  class basic_teebuf
        : public std::basic_streambuf<charT, traits>
    {
        std::basic_streambuf<charT, traits>* sb1_;
        std::basic_streambuf<charT, traits>* sb2_;

        int overflow(int c) {
            bool rc(true);
            if (!traits::eq_int_type(traits::eof(), c)) {
                traits::eq_int_type(this->sb1_->sputc(c), traits::eof())
                    && (rc = false);
                traits::eq_int_type(this->sb2_->sputc(c), traits::eof())
                    && (rc = false);
            }
            return rc ? traits::not_eof(c) : traits::eof();
        }
        int sync() {
            bool rc(false);
            this->sb1_->pubsync() != -1 || (rc = false);
            this->sb2_->pubsync() != -1 || (rc = false);
            return rc ? -1 : 0;
        }
    public:
        basic_teebuf(std::basic_streambuf<charT, traits>* sb1, std::basic_streambuf<charT, traits>* sb2)
            : sb1_(sb1), sb2_(sb2) {
        }
        
    };

    typedef basic_teebuf<char> teebuf;

    template <class charT, class traits = std::char_traits<charT>> class basic_oteestream
        : private virtual basic_teebuf<charT, traits>
        , public std::basic_ostream<charT, traits> {
    public:
        basic_oteestream(std::basic_ostream<charT, traits>& out1, std::basic_ostream<charT, traits>& out2)
            : basic_teebuf<charT, traits>(out1.rdbuf(), out2.rdbuf())
            , std::basic_ostream<charT, traits>(this) {
            this->init(this);
        }
    };

    typedef basic_oteestream<char> oteestream;

    teebuf *tee_cout(std::filesystem::path& opath);
    teebuf *tee_cerr(std::filesystem::path& opath);
    teebuf *tee_cout_cerr(std::filesystem::path& opath);
}
