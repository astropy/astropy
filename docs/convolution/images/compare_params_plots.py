# This code is meant to be run interactively (%run -i equivalent) in an
# environment where compare_params.py has already been run.  It has been split
# out to reduce the amount of "do-nothing" code shown in the docs.
#
# Now we do a bunch of plots.  In the first two plots, the originally masked
# values are marked with red X's
plt.figure(1, figsize=(12,12)).clf()
ax1 = plt.subplot(2,2,1)
im = ax1.imshow(img, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
y,x = np.where(np.isnan(img))
ax1.set_autoscale_on(False)
ax1.plot(x, y, 'rx', markersize=4)
ax1.set_title("Original")
ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax2 = plt.subplot(2,2,2)
im = ax2.imshow(scipy_conv, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax2.set_autoscale_on(False)
ax2.plot(x, y, 'rx', markersize=4)
ax2.set_title("Scipy")
ax2.set_xticklabels([])
ax2.set_yticklabels([])

ax3 = plt.subplot(2,2,3)
im = ax3.imshow(scipy_conv_zerod, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax3.set_title("Scipy nan->zero")
ax3.set_xticklabels([])
ax3.set_yticklabels([])


ax4 = plt.subplot(2,2,4)
im = ax4.imshow(astropy_conv, vmin=-2., vmax=2.e1, origin='lower',
                    interpolation='nearest', cmap='viridis')
ax4.set_title("Default astropy")
ax4.set_xticklabels([])
ax4.set_yticklabels([])


# we make a second plot of the amplitudes vs offset position to more clearly
# illustrate the value differences
plt.figure(2).clf()
plt.plot(img[:,25], label='input', drawstyle='steps-mid', linewidth=2, alpha=0.5)
plt.plot(scipy_conv[:,25], label='scipy', drawstyle='steps-mid',
         linewidth=2, alpha=0.5, marker='s')
plt.plot(scipy_conv_zerod[:,25], label='scipy nan->zero', drawstyle='steps-mid',
         linewidth=2, alpha=0.5, marker='s')
plt.plot(astropy_conv[:,25], label='astropy', drawstyle='steps-mid', linewidth=2, alpha=0.5)
plt.ylabel("Amplitude")
plt.ylabel("Position Offset")
plt.legend(loc='best')
plt.show()
