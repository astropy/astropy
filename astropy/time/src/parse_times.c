#include <stdio.h>
#include <string.h>

// ASCII codes for '0' and '9'
const char char_zero = 48;
const char char_nine = 57;

// Distutils on Windows automatically exports ``PyInit__parse_times``,
// create dummy to prevent linker complaining about missing symbol.
// Based on convolution/src/convolve.c.
#if defined(_MSC_VER)
void PyInit__parse_times(void)
{
    return;
}
#endif

int parse_int_from_char_array(char *chars, int str_len,
                              char delim, int idx0, int idx1,
                              int *val)
// Parse integer from positions idx0:idx1 (inclusive) within chars, optionally
// starting with a delimiter.
//
// Example: "2020-01-24"
//                  ^^^
//           0123456789
//
// int day, status;
// status = parse_int_from_char_array("2020-01-24", &day, 10, '-', 7, 9);
//
// Inputs:
//  char *chars: time string
//  int str_len: length of *chars string
//  char delim: optional character at position idx0 when delim > 0
//  int idx0: start index for parsing integer
//  int idx1: stop index (inclusive) for parsing integer
//
// Output:
//  int *val: output value
//
// Returns:
//  int status:
//    0: OK
//    1: String ends at the beginning of requested value
//    2: String ends in the middle of requested value
//    3: Required delimiter character not found
//    4: Non-digit found where digit (0-9) required
{
    int mult = 1;
    char digit;
    char ch;
    int ii;

    // Check if string ends (has 0x00) before str_len. Require that this segment
    // of the string is entirely contained in the string (idx1 < str_len),
    // remembering that idx1 is inclusive and counts from 0.
    if (idx1 < str_len) {
        for (ii = idx0; ii <= idx1; ii++) {
            if (chars[ii] == 0) {
                str_len = ii;
                break;
            }
        }
    }
    // String ends before the beginning of requested value,
    // e.g. "2000-01" (str_len=7) for day (idx0=7). This is OK in some
    // cases, e.g. before hour (2000-01-01).
    if (idx0 >= str_len) {
        return 1;
    }

    // String ends in the middle of requested value. This implies a badly
    // formatted time.
    if (idx1 >= str_len) {
        return 2;
    }

    // Look for optional delimiter character, e.g. ':' before minute. If delim == 0
    // then no character is required.
    if (delim > 0) {
        // Required start character not found.
        if (chars[idx0] != delim) {
            return 3;
        }
        idx0 += 1;
    }

    // Build up the value using reversed digits
    *val = 0;
    for (ii = idx1; ii >= idx0; ii--)
    {
        ch = chars[ii];
        if (ch < char_zero || ch > char_nine) {
            // Not a digit, implying badly formatted time.
            return 4;
        }
        digit = ch - char_zero;
        *val += digit * mult;
        mult *= 10;
    }

    return 0;
}

int parse_frac_from_char_array(char *chars, int str_len, char delim, int idx0,
                               double *val)
// Parse trailing fraction starting from position idx0 in chars.
//
// Example: "2020-01-24T12:13:14.5556"
//                              ^^^^^
//           012345678901234567890123
//
// int status;
// float frac;
// status = parse_frac_from_char_array("2020-01-24T12:13:14.5556", &frac, 24, '.', 19);
//
// Inputs:
//  char *chars: time string
//  int str_len: length of *chars string
//  char delim: optional character at position idx0 when delim > 0
//  int idx0: start index for parsing integer
//
// Output:
//  double *val: output value
//
// Returns:
//  int status:
//    0: OK
//    1: String ends at the beginning of requested value
//    3: Required delimiter character not found
//    4: Non-digit found where digit (0-9) required
{
    double mult = 0.1;
    char digit;
    char ch;
    int ii;

    *val = 0.0;

    // String ends at exactly before the beginning of requested fraction.
    // e.g. "2000-01-01 12:13:14". Fraction value is zero.
    if (idx0 == str_len) {
        return 1;
    }

    // Look for optional delimiter character, e.g. '.' before fraction. If delim == 0
    // then no character is required. This can happen for unusual formats like
    // Chandra GRETA time yyyyddd.hhmmssfff.
    if (delim > 0) {
        // Required start character not found.
        if (chars[idx0] != delim) {
            return 3;
        }
        idx0 += 1;
    }

    for (ii = idx0; ii < str_len; ii++)
    {
        ch = chars[ii];
        if (ch < char_zero || ch > char_nine) {
            // Not a digit, implying badly formatted time.
            return 4;
        }
        digit = ch - char_zero;
        *val += digit * mult;
        mult /= 10.0;
    }
    return 0;
}

static inline int is_leap_year (int year)
// Determine if year is a leap year.
// Inspired by from https://stackoverflow.com/questions/17634282
{
  return ((year & 3) == 0)
          && ((year % 100 != 0)
              || (((year / 100) & 3) == 0));
}

int convert_day_of_year_to_month_day(int year, int day_of_year, int *month, int *day_of_month)
// Convert year and day_of_year into month, day_of_month
// Inspired by from https://stackoverflow.com/questions/17634282, determine
{
    int leap_year = is_leap_year(year) ? 1 : 0;
    int days_in_year = leap_year ? 366 : 365;
    const unsigned short int _mon_yday_normal[13] =
        { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
    const unsigned short int _mon_yday_leap[13] =
        { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };
    const unsigned short int *mon_yday = leap_year ? _mon_yday_leap :_mon_yday_normal;
    int mon;

    if (day_of_year < 1 || day_of_year > days_in_year) {
        // Error in day_of_year
        return 5;
    }

    for (mon = 1; mon <= 12; mon++) {
        if (day_of_year <= mon_yday[mon]) {
            *month = mon;
            *day_of_month = day_of_year - mon_yday[mon - 1];
            break;
        }
    }

    return 0;
}

int parse_ymdhms_times(char *times, int n_times, int max_str_len, int has_day_of_year,
                   char *delims, int *starts, int *stops, int *break_allowed,
                   int *years, int *months, int *days, int *hours,
                   int *minutes, double *seconds)
// Parse a string time in `chars` which has year, (month, day | day_of_year),
// hour, minute, seconds components.
//
// Examples: "2020-01-24T12:13:14.5556", "2020:123:12:13:14.5556"
//
// Inputs:
//  char *times: time characters (flattened n_times x max_str_len array)
//  int n_times: number of time strings (each max_str_len long)
//  int max_str_len: max length of string (may be null-terminated before this)
//  int has_day_of_year: time includes day-of-year instead of month, day-of-month
//  char *delims: array of delimiters preceding yr, mon, day, hr, min, isec, frac
//      components. Value of 0 means no preceding delimiter.
//  int *starts, *stop: arrays of start/stop indexes into time string.
//  int *break_allowed: if true (1) then the time string can legally end just
//      before the corresponding component (e.g. "2000-01-01" is a valid time but
//      "2000-01-01 12" is not).
//
// Outputs:
//  int *year, *month, *day, *hour, *minute: output components (n_times long)
//  double *second: output seconds (n_times long)
//
// Returns:
//  int status:
//    0: OK
//    1: String ends at the beginning of requested value
//    2: String ends in the middle of requested value
//    3: Required delimiter character not found
//    4: Non-digit found where digit (0-9) required
//    5: Bad day of year
{
    int str_len;
    int status = 0;
    int isec;
    double frac;
    char *time;
    int *year, *month, *day, *hour, *minute;
    double *second;
    int i, ii;

    for (ii = 0; ii < n_times; ii++)
    {
        time = times + ii * max_str_len;
        year = years + ii;
        month = months + ii;
        day = days + ii;
        hour = hours + ii;
        minute = minutes + ii;
        second = seconds + ii;

        // Initialize default values
        *month = 1;
        *day = 1;
        *hour = 0;
        *minute = 0;
        *second = 0.0;

        // Parse "2000-01-12 13:14:15.678"
        //        01234567890123456789012

        // Check for null termination before max_str_len. If called using a contiguous
        // numpy 2-d array of chars there may or may not be null terminations.
        str_len = max_str_len;
        for (i = 0; i < max_str_len; i++) {
            if (time[i] == 0) {
                str_len = i;
                break;
            }
        }

        // Get each time component year, month, day, hour, minute, isec, frac
        status = parse_int_from_char_array(time, str_len, delims[0], starts[0], stops[0], year);
        if (status) {
            if (status == 1 && break_allowed[0]) { continue; }
            else { return status; }
        }

        // Optionally parse month
        if (! has_day_of_year) {
            status = parse_int_from_char_array(time, str_len, delims[1], starts[1], stops[1], month);
            if (status) {
                if (status == 1 && break_allowed[1]) { continue; }
                else { return status; }
            }
        }

        // This might be day-of-month or day-of-year
        status = parse_int_from_char_array(time, str_len, delims[2], starts[2], stops[2], day);
        if (status) {
            if (status == 1 && break_allowed[2]) { continue; }
            else { return status; }
        }

        if (has_day_of_year) {
            // day contains day of year at this point, but convert it to day of month
            status = convert_day_of_year_to_month_day(*year, *day, month, day);
            if (status) {
                return status;
            }
        }

        status = parse_int_from_char_array(time, str_len, delims[3], starts[3], stops[3], hour);
        if (status) {
            if (status == 1 && break_allowed[3]) { continue; }
            else { return status; }
        }

        status = parse_int_from_char_array(time, str_len, delims[4], starts[4], stops[4], minute);
        if (status) {
            if (status == 1 && break_allowed[4]) { continue; }
            else { return status; }
        }

        status = parse_int_from_char_array(time, str_len, delims[5], starts[5], stops[5], &isec);
        if (status) {
            if (status == 1 && break_allowed[5]) { continue; }
            else { return status; }
        }

        status = parse_frac_from_char_array(time, str_len, delims[6], starts[6], &frac);
        if (status) {
            if (status != 1 || ! break_allowed[6]) { return status; }
        }

        *second = isec + frac;
    }

    return 0;
}

int check_unicode(char *chars, int n_unicode_char)
// Check if *chars is pure ASCII, assuming input is UTF-32
{
    char *ch;
    int ii;

    ch = chars;
    for (ii = 0; ii < n_unicode_char; ii++)
    {
        ch++;
        if (*ch++) return 1;
        if (*ch++) return 1;
        if (*ch++) return 1;
    }
    return 0;

}
