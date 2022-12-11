#include "wirth_select.h"
#include <math.h>
#include <stdlib.h>

void compute_sigma_clipped_bounds(double data_buffer[], int count, int use_median,
                                  int use_mad_std, int maxiters, double sigma_lower,
                                  double sigma_upper, double *lower_bound,
                                  double *upper_bound, double  mad_buffer[]) {

  double mean, std, median = 0, cen;
  int i, new_count, iteration = 0;

  while (1) {

    if (use_median || use_mad_std) {
      median = wirth_median(data_buffer, count);
    }

    // Note that we don't use an else clause here, because the mean
    // might be needed even if use_median is used, but use_mad_std is not.
    if (!use_median || !use_mad_std) {
      mean = 0;
      for (i = 0; i < count; i++) {
        mean += data_buffer[i];
      }
      mean /= count;
    }

    if (use_median) {
      cen = median;
    } else {
      cen = mean;
    }

    if (use_mad_std) {

      for (i = 0; i < count; i++) {
        mad_buffer[i] = fabs(data_buffer[i] - median);
      }
      std = wirth_median(mad_buffer, count) * 1.482602218505602;

    } else {

      std = 0;
      for (i = 0; i < count; i++) {
        std += pow(mean - data_buffer[i], 2);
      }
      std = sqrt(std / count);

    }

    *lower_bound = cen - sigma_lower * std;
    *upper_bound = cen + sigma_upper * std;

    // We now exclude values from the buffer using these
    // limits and shift values so that we end up with a
    // packed array of 'valid' values
    new_count = 0;
    for (i = 0; i < count; i++) {
      if (data_buffer[i] >= *lower_bound && data_buffer[i] <= *upper_bound) {
        data_buffer[new_count] = data_buffer[i];
        new_count += 1;
      }
    }

    if (new_count == count)
      return;

    count = new_count;

    iteration += 1;

    if (maxiters != -1 && iteration >= maxiters)
      return;
  }
}
