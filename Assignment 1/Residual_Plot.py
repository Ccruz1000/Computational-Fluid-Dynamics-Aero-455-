import matplotlib.pyplot as plt
import os
num_panels = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
residual1 = [1.1852562233705655e-15, 1.3322676295501878e-15, 1.2212453270876722e-15, 3.885780586188048e-16,
             1.609823385706477e-15, 1.2212453270876722e-15, 2.546574062733953e-15, 3.969047313034935e-15,
             -4.0505793164058446e-15, 2.424276057677588e-15, -1.0148132334464322e-16, 2.424276057677588e-16,
             1.2445556737961105e-15]
plt.figure(3)
plt.plot(num_panels, residual1, color='r', label='Residual')
plt.title('Source Panel Method Residual')
plt.xlabel('Number of Panels')
plt.ylabel('Residual')
plt.legend(loc='best')
current_path = os.getcwd()
plot_folder = current_path + '/' + 'Residual'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder, exist_ok=True)
plt.savefig(plot_folder + '/Residual.png', bbox_extra_artists='legend_outside')
plt.show()
# plt.close()