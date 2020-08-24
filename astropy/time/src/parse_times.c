#include <stdio.h>
#include <string.h>

const char char_zero = 48;
const char char_nine = 57;

int parse_int_from_char_array(char *chars, int *val, int str_len,
                              char char_start, int idx0, int idx1)
// Parse an integer from positions idx0:idx1 (inclusive) within chars.
//
// Example: "2020-01-24"
//                  ^^^
//           0123456789
//
// int day, status;
// status = parse_int_from_char_array("2020-01-24", &day, 10, '-', 7, 9);
//
// Args:
//  char *chars: time string
//  int *val: output value
//  int str_len: length of *chars string
//  char char_start: optional character at position idx0 when char_start > 0
//  int idx0: start index for parsing integer
//  int idx1: stop index (inclusive) for parsing integer
{
    int mult = 1;
    char digit;
    char ch;
    int status = 0;

    // Check if string ends (has 0x00) before str_len. Require that this segment
    // of the string is entirely contained in the string (idx1 < str_len),
    // remembering that idx1 is inclusive and counts from 0.
    if (idx1 < str_len) {
        for (size_t i = idx0; i <= idx1; i++) {
            if (chars[i] == 0) {
                str_len = i;
                break;
            }
        }
    }
    // String ends before the beginning of requested value,
    // e.g. "2000-01" (str_len=7) for day (idx0=7). This is OK in some
    // cases, e.g. before hour (2000-01-01).
    if (idx0 >= str_len) {
        return -1;
    }

    // String ends in the middle of requested value. This implies a badly
    // formatted time.
    if (idx1 >= str_len) {
        return -2;
    }

    // Look for optional start character, e.g. ':' before minute. If char_start == 0
    // then no character is required.
    if (char_start > 0) {
        // Required start character not found.
        if (chars[idx0] != char_start) {
            return -3;
        }
        idx0 += 1;
    }

    // Build up the value using reversed digits
    *val = 0;
    for (int ii = idx1; ii >= idx0; ii--)
    {
        ch = chars[ii];
        if (ch < char_zero || ch > char_nine) {
            // Not a digit, implying badly formatted time.
            return -4;
        }
        digit = ch - char_zero;
        *val += digit * mult;
        mult *= 10;
    }

    return 0;
}

int parse_frac_from_char_array(char *chars, double *val,
                               int str_len, char char_start, int idx0)
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
// Args:
//  char *chars: time string
//  double *val: output fraction value
//  int str_len: length of *chars string
//  char char_start: optional character at position idx0 when char_start > 0
//  int idx0: start index for parsing fraction
{
    double mult = 0.1;
    char digit;
    char ch;
    int status = 0;

    *val = 0.0;

    // String ends at exactly before the beginning of requested fraction.
    // e.g. "2000-01-01 12:13:14". Fraction value is zero.
    if (idx0 == str_len) {
        return 0;
    }

    // Look for optional start character, e.g. '.' before fraction. If char_start == 0
    // then no character is required. This can happen for unusual formats like
    // Chandra GRETA time yyyyddd.hhmmssfff.
    if (char_start > 0) {
        // Required start character not found.
        if (chars[idx0] != char_start) {
            return -3;
        }
        idx0 += 1;
    }

    for (size_t ii = idx0; ii < str_len; ii++)
    {
        ch = chars[ii];
        if (ch < char_zero || ch > char_nine) {
            // Not a digit, implying badly formatted time.
            return -4;
        }
        digit = ch - char_zero;
        *val += digit * mult;
        mult /= 10.0;
    }
    return 0;
}

int parse_ymdhms_times(char *times, int n_times, int max_str_len,
                   char *delims, int *starts, int *stops, int *break_allowed,
                   int *years, int *months, int *days, int *hours,
                   int *minutes, double *seconds)
// Parse an ISO time in `chars`.
//
// Example: "2020-01-24T12:13:14.5556"
//
// Args:
//  char *times: time characters (flattened n_times x max_str_len array)
//  int n_times: number of time strings (each max_str_len long)
//  int max_str_len: max length of string (may be null-terminated before this)
//  int *year, *month, *day, *hour, *minute: output components (n_times long)
//  double *second: output seconds (n_times long)
//
// Returns:
//  int status: 0 for OK, < 0 for not OK
{
    int str_len;
    int status = 0;
    int isec;
    double frac;
    char sep = ' ';
    char *time;
    int *year, *month, *day, *hour, *minute;
    double *second;

    for (size_t ii = 0; ii < n_times; ii++)
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
        for (size_t i = 0; i < max_str_len; i++) {
            if (time[i] == 0) {
                str_len = i;
                break;
            }
        }

        status = parse_int_from_char_array(time, year, str_len, delims[0], starts[0], stops[0]);
        if (status < 0) {
            if (status == -1 && break_allowed[0]) { continue; }
            else { return status; }
        }

        status = parse_int_from_char_array(time, month, str_len, delims[1], starts[1], stops[1]);
        if (status < 0) {
            if (status == -1 && break_allowed[1]) { continue; }
            else { return status; }
        }

        status = parse_int_from_char_array(time, day, str_len, delims[2], starts[2], stops[2]);
        if (status < 0) {
            if (status == -1 && break_allowed[2]) { continue; }
            else { return status; }
        }

        status = parse_int_from_char_array(time, hour, str_len, delims[3], starts[3], stops[3]);
        if (status < 0) {
            if (status == -1 && break_allowed[3]) { continue; }
            else { return status; }
        }

        status = parse_int_from_char_array(time, minute, str_len, delims[4], starts[4], stops[4]);
        if (status < 0) {
            if (status == -1 && break_allowed[4]) { continue; }
            else { return status; }
        }

        status = parse_int_from_char_array(time, &isec, str_len, delims[5], starts[5], stops[5]);
        if (status < 0) {
            if (status == -1 && break_allowed[5]) { continue; }
            else { return status; }
        }

        status = parse_frac_from_char_array(time, &frac, str_len, delims[6], starts[6]);
        if (status < 0) { return status; }

        *second = isec + frac;
    }

    return 0;
}

int check_unicode(char *chars, int n_unicode_char) {
    char *ch;

    ch = chars;
    for (size_t i = 0; i < n_unicode_char; i++)
    {
        ch++;
        if (*ch++) return -1;
        if (*ch++) return -1;
        if (*ch++) return -1;
    }
    return 0;

}
