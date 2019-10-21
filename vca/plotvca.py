##########################
#Plot the VCA by channel

specin=vca_array[:,0]
thickness=vca_array[:,1]
upper_err=vca_array[:,0]+vca_array[:,2]
lower_err=vca_array[:,0]-vca_array[:,2]

plt.title('VCA:' + name)
plt.xlabel('Channel Thickness (km/s)')
plt.ylabel('Spectral Index')
plt.tight_layout()
plt.fill_between(thickness, upper_err, lower_err)
plt.scatter(thickness,specin)
plt.tight_layout()
plt.savefig(figsaveloc + '.png', format='png')
