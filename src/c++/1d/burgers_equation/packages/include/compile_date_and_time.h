//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D Burgers equation, fully implicit with delta-formulation and Modified Newton iteration 
//    Copyright (C) 2026 Jan Mooiman
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------------

#ifndef __COMPILE_DATE_AND_TIME__
#define __COMPILE_DATE_AND_TIME__

constexpr std::string_view getMonthNumber(std::string_view month) {
    return month == "Jan" ? "01" :
           month == "Feb" ? "02" :
           month == "Mar" ? "03" :
           month == "Apr" ? "04" :
           month == "May" ? "05" :
           month == "Jun" ? "06" :
           month == "Jul" ? "07" :
           month == "Aug" ? "08" :
           month == "Sep" ? "09" :
           month == "Oct" ? "10" :
           month == "Nov" ? "11" :
           month == "Dec" ? "12" :
                            "??";
}

constexpr char toDigit(char c) {
    return (c >= '0' && c <= '9') ? c : ' ';
}

constexpr std::string_view compileYear() {
    constexpr std::string_view date = __DATE__;  // "Jul 11 2025"

    constexpr std::string_view monthStr = date.substr(0, 3);
    constexpr std::string_view dayStr = (date[4] == ' ') ? date.substr(5, 1) : date.substr(4, 2);
    constexpr std::string_view yearStr = date.substr(7, 4);

    constexpr std::string_view monthNum = getMonthNumber(monthStr);

    // NOTE: Cannot construct std::string_view concatenation directly in constexpr yet
    return yearStr; // For now just returning year as proof of concept
}
constexpr std::string_view compileMonth() {
    constexpr std::string_view date = __DATE__;  // "Jul 11 2025"

    constexpr std::string_view monthStr = date.substr(0, 3);
    constexpr std::string_view dayStr = (date[4] == ' ') ? date.substr(5, 1) : date.substr(4, 2);
    constexpr std::string_view yearStr = date.substr(7, 4);

    constexpr std::string_view monthNum = getMonthNumber(monthStr);

    // NOTE: Cannot construct std::string_view concatenation directly in constexpr yet
    return monthNum; // For now just returning year as proof of concept
}
constexpr std::string_view compileDay() {
    constexpr std::string_view date = __DATE__;  // "Jul 11 2025"

    constexpr std::string_view monthStr = date.substr(0, 3);
    constexpr std::string_view dayStr = (date[4] == ' ') ? date.substr(5, 1) : date.substr(4, 2);
    constexpr std::string_view yearStr = date.substr(7, 4);

    constexpr std::string_view monthNum = getMonthNumber(monthStr);

    // NOTE: Cannot construct std::string_view concatenation directly in constexpr yet
    return dayStr; // For now just returning year as proof of concept
}
std::string compileDateTime()
{
    std::string str1(compileYear());
    std::string str2(compileMonth());
    std::string str3(compileDay());
    if (str3.size() == 1)
    {
        str3.resize(2);
        str3[1] = str3[0];
        str3[0] = '0';
    }
    return str1 + "-" + str2 + "-" + str3 + " " + __TIME__;
}
#endif __COMPILE_DATE_AND_TIME__
