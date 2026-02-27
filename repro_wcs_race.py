
import threading
import numpy as np
import time
from astropy.wcs import WCS
from astropy.io import fits
import warnings

# Suppress warnings to keep output clean
warnings.filterwarnings("ignore")

def make_wcs():
    # Create a simple WCS with some distortion to ensure we hit the C code paths
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [100.0, 45.0]
    w.wcs.crpix = [500.0, 500.0]
    w.wcs.cdelt = [-0.01, 0.01]
    
    # Add some SIP distortion to make it interesting (and use more C code)
    # This requires creating a header-based WCS usually, but we can try setting sip directly if easier.
    # For now, a basic WCS is often enough to trigger weirdness if shared.
    w.wcs.set()
    return w

def worker(w, thread_id, n_iter):
    # Create random coordinates
    coords = np.random.rand(100, 2) * 1000
    
    for i in range(n_iter):
        try:
            # Race condition: accessing the same w.wcs struct
            # wcs_pix2world writes to w.wcs internal buffers
            res = w.all_pix2world(coords, 0)
            
            # Sanity check (optional, mostly looking for segfaults)
            if not np.isfinite(res).all():
                print(f"Thread {thread_id}: NaN result!")
        except Exception as e:
            print(f"Thread {thread_id} error: {e}")

def run_race():
    print("Initializing WCS...")
    try:
        w = make_wcs()
    except Exception as e:
        print(f"Failed to create WCS: {e}")
        return

    n_threads = 16
    n_iter = 1000
    threads = []
    
    print(f"Starting {n_threads} threads for {n_iter} iterations...")
    
    for i in range(n_threads):
        t = threading.Thread(target=worker, args=(w, i, n_iter))
        threads.append(t)
        t.start()
        
    for t in threads:
        t.join()
        
    print("Finished without segfault (lucky?).")

if __name__ == "__main__":
    run_race()
