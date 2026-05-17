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
#endif __COMPILE_DATE_AND_TIME__
