//
// Programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//
//
//    Solving the 2D shallow water equations, fully implicit with delta-formuation and Modified Newton iteration 
//    Copyright (C) 2025 Jan Mooiman
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
#ifndef __PERF_TIMER_H__
#define __PERF_TIMER_H__

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <iterator>

#define DO_TIMING 1
#if DO_TIMING == 1
#   define TOSTRING(x) #x
#   define START_TIMER(x) tt = TestTimer::get_instance(); tt->start(TOSTRING(x), __FILE__, __FUNCTION__, __LINE__);
#   define START_TIMERN(x) TestTimer * tt = TestTimer::get_instance(); tt->start(TOSTRING(x), __FILE__, __FUNCTION__, __LINE__);
#   define CLEAR_TIMER() if (tt->get_count() != 0){ tt->clear(); }
#   define DELETE_TIMER() delete tt ;
#   define DELETE_TIMERN() TestTimer * tt = TestTimer::get_instance(); delete tt ;
#   define STOP_TIMER(x) tt->stop(TOSTRING(x));
#   define PRINT_TIMER(x) tt->log(x);  // parameter x is a filename
#   define PRINT_TIMERN(x) TestTimer * tt = TestTimer::get_instance(); tt->log(x);  // paramter x is a filename
#else
#   define START_TIMERN(x)
#   define START_TIMER(x)
#   define CLEAR_TIMER()
#   define DELETE_TIMER()
#   define DELETE_TIMERN()
#   define STOP_TIMER(x)
#   define PRINT_TIMER(x)
#   define PRINT_TIMERN(x)
#   define PRINT_TIMERS(x)
#endif

class Timing {
public:
    Timing() = default;
    Timing(std::string name,
        std::string fileName,
        std::string functionName,
        int lineNumber) :
        level(level), name(std::move(name)), filename(std::move(fileName)), functionname(std::move(functionName)), line_number(lineNumber)
    {
        start = std::chrono::high_resolution_clock::now();
        elapse_time = 0.0;
        level = 0;
        nr_calls = 1;
        running = false;
    }

    bool running;
    std::string name;
    int level;
    int line_number;
    int nr_calls;
    double elapse_time;
    std::string filename;
    std::string functionname;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point stop;
};



class TestTimer
{
public:
    static TestTimer* get_instance()
    {
        if (obj == nullptr)
        {
            obj = new TestTimer(std::string("root"),
                std::string("---"),
                std::string("---"),
                0);
            obj->m_level = 0;
        }
        object_count = 1;
        return obj;
    }

    int get_count();
    void log(std::string);
    std::vector<Timing> timing;

    TestTimer(const std::string& name, 
        const std::string& fileName,
        const std::string& functionName,
        int lineNumber);

    ~TestTimer() = default;

    void clear();
    void start(const std::string& name, const std::string& fileName, const std::string& functionName, int lineNumber);
    void stop(const std::string& name);

    private:
        static TestTimer * obj;
        static int object_count;

        int m_level;
        std::string m_name;
        std::string m_filename;
        std::string m_functionname;
        int m_line_number;
};
#endif  // __PERF_TIMER_H__
