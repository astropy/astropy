
/*
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 */

#define ELEM_SWAP(a, b)                                                        \
  {                                                                            \
    register double t = (a);                                                   \
    (a) = (b);                                                                 \
    (b) = t;                                                                   \
  }

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median.

                Reference:

                  Author: Wirth, Niklaus
                   Title: Algorithms + data structures = programs
               Publisher: Englewood Cliffs: Prentice-Hall, 1976
    Physical description: 366 p.
                  Series: Prentice-Hall Series in Automatic Computation

 ---------------------------------------------------------------------------*/

double kth_smallest(double a[], int n, int k) {
  register int i, j, l, m;
  register double x;

  l = 0;
  m = n - 1;
  while (l < m) {
    x = a[k];
    i = l;
    j = m;
    do {
      while (a[i] < x)
        i++;
      while (x < a[j])
        j--;
      if (i <= j) {
        ELEM_SWAP(a[i], a[j]);
        i++;
        j--;
      }
    } while (i <= j);
    if (j < k)
      l = i;
    if (k < i)
      m = j;
  }
  return a[k];
}

double wirth_median(double a[], int n) {
  if (n % 2 == 0) {
    return 0.5 * (kth_smallest(a, n, n / 2) + kth_smallest(a, n, n / 2 - 1));
  } else {
    return kth_smallest(a, n, (n - 1) / 2);
  }
}
