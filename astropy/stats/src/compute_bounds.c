#include "wirth_select.h"
#include <math.h>

void compute_sigma_clipped_bounds(double buffer[], int count, int use_median,
                                  int maxiters, double sigma_lower,
                                  double sigma_upper, double *lower_bound,
                                  double *upper_bound) {

  double mean, std, median;
  int i, new_count, iteration = 0;

  while (1) {

    mean = 0;
    for (i = 0; i < count; i++) {
      mean += buffer[i];
    }
    mean /= count;

    std = 0;
    for (i = 0; i < count; i++) {
      std += pow(mean - buffer[i], 2);
    }
    std = sqrt(std / count);

    if (use_median) {

      median = wirth_median(buffer, count);

      *lower_bound = median - sigma_lower * std;
      *upper_bound = median + sigma_upper * std;

    } else {

      *lower_bound = mean - sigma_lower * std;
      *upper_bound = mean + sigma_upper * std;
    }

    // We now exclude values from the buffer using these
    // limits and shift values so that we end up with a
    // packed array of 'valid' values
    new_count = 0;
    for (i = 0; i < count; i++) {
      if (buffer[i] >= *lower_bound && buffer[i] <= *upper_bound) {
        buffer[new_count] = buffer[i];
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
