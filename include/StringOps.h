//=================================
// include guard
#ifndef __STRINGOPS_INCLUDED__
#define __STRINGOPS_INCLUDED__

// #include <algorithm> 
// #include <cctype>
// #include <locale>

// // trim from start (in place)
// static inline void ltrim(std::string &s) {
//     s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
//         return !std::isspace(ch);
//     }));
// }

// // trim from end (in place)
// static inline void rtrim(std::string &s) {
//     s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
//         return !std::isspace(ch);
//     }).base(), s.end());
// }

// // trim from both ends (in place)
// static inline void trim(std::string &s) {
//     ltrim(s);
//     rtrim(s);
// }

// // trim from start (copying)
// static inline std::string ltrim_copy(std::string s) {
//     ltrim(s);
//     return s;
// }

// // trim from end (copying)
// static inline std::string rtrim_copy(std::string s) {
//     rtrim(s);
//     return s;
// }

// // trim from both ends (copying)
// static inline std::string trim_copy(std::string s) {
//     trim(s);
//     return s;
// }


#include <string>

// modifies input string, returns input

inline std::string& trim_left_in_place(std::string& str) {
    size_t i = 0;
    while(i < str.size() && isspace(str[i])) { ++i; };
    return str.erase(0, i);
}

inline std::string& trim_right_in_place(std::string& str) {
    size_t i = str.size();
    while(i > 0 && isspace(str[i - 1])) { --i; };
    return str.erase(i, str.size());
}

inline std::string& trim_in_place(std::string& str) {
    return trim_left_in_place(trim_right_in_place(str));
}

// returns newly created strings

inline std::string trim_right(std::string str) {
    return trim_right_in_place(str);
}

inline std::string trim_left(std::string str) {
    return trim_left_in_place(str);
}

inline std::string trim(std::string str) {
    return trim_left_in_place(trim_right_in_place(str));
}

#endif