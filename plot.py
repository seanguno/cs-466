import matplotlib.pyplot as plt
import numpy as np

# plot time analysis
# x-axis = average length of the two aligned sequences
# y-axis = average time of alignment (alignment performed 5 times, average time shown)
nw_time = [(19.0, 0.00035889982245862484), (479.0, 0.2004694000352174), (1001.5, 0.9101809000130743), (1763.5, 2.6902367000002414), (2229.0, 3.9246781999245286), (6403.5, 36.608175799949095), (9307.0, 75.01152989990078), (16569.0, 256.68853710009716)]
hb_time = [(19.0, 0.00038759992457926273), (479.0, 0.14335140003822744), (1001.5, 0.6156985000707209), (1763.5, 1.902610900113359), (2229.0, 3.091165600111708), (6403.5, 26.724719499936327), (9307.0, 58.57902719988488), (16569.0, 192.69252449995838)] 

nw_time_extended = [(19.0, 0.0004990999586880207), (479.0, 0.27225330006331205), (1001.5, 1.2385339997708797), (1763.5, 3.7759147002361715), (2229.0, 4.180357899982482), (6403.5, 34.75250760000199), (9307.0, 74.01478530000895), (16569.0, 236.11566699994728), (1000, 0.8331424002535641), (2000, 3.3964247000403702), (3000, 7.649441899731755), (4000, 13.700268200132996), (5000, 21.75459419982508), (6000, 32.62650189967826), (7000, 41.84669720008969), (8000, 55.991296799853444), (9000, 72.49378219991922), (10000, 86.82874789973721), (11000, 102.13932760013267), (12000, 122.52714310027659), (13000, 143.80858709989116), (14000, 219.46323610004038), (15000, 275.8117311000824)]
hb_time_extended = [(19.0, 0.0006675003096461296), (479.0, 0.21720070019364357), (1001.5, 0.9112697001546621), (1763.5, 2.3741548000834882), (2229.0, 3.608909399714321), (6403.5, 36.08737460011616), (9307.0, 66.23498589964584), (16569.0, 223.51840800000355), (1000, 0.7073352998122573), (2000, 2.8966955998912454), (3000, 6.6584540996700525), (4000, 12.199166500009596), (5000, 19.02023210003972), (6000, 27.1168828997761), (7000, 37.55801279982552), (8000, 50.481557599734515), (9000, 65.34488129988313), (10000, 80.00038050021976), (11000, 95.23744589975104), (12000, 113.36187560018152), (13000, 135.7678286000155), (14000, 187.17405309993774), (15000, 186.1493144002743)]
x, y = zip(*nw_time_extended)
x2, y2 = zip(*hb_time_extended)

# fit quadratic polynomial to both
coeffs_nw = np.polyfit(x, y, 2)
poly_nw = np.poly1d(coeffs_nw)
print(f'NW runtime quadratic fit coefficients: {coeffs_nw}')

coeffs_hb = np.polyfit(x2, y2, 2)
poly_hb = np.poly1d(coeffs_hb)
print(f'Hirschberg runtime quadratic fit coefficients: {coeffs_hb}')

# x-vals for polynomial plotting
x_fit = np.linspace(min(min(x), min(x2)), max(max(x), max(x2)), 100)

# plotting
plt.figure(1)
plt.scatter(x, y, color='blue', label='Needleman-Wunsch')
plt.scatter(x2, y2, color='red', label='Hirschberg')

plt.plot(x_fit, poly_nw(x_fit), color='blue', linestyle='--', label = 'Quadratic Fit for Needleman-Wunsch')
plt.plot(x_fit, poly_hb(x_fit), color='red', linestyle='--', label = 'Quadratic Fit for Hirschberg')

plt.xlabel('Average Sequence Length')
plt.ylabel('Average Alignment Time (sec)')
plt.legend()

# calculating R^2 values
y_pred = poly_nw(x)
y2_pred = poly_hb(x2)
SS_res = np.sum((y - y_pred) ** 2)
SS2_res = np.sum((y2 - y2_pred) ** 2)
SS_total = np.sum((y - np.mean(y)) ** 2)
SS2_total = np.sum((y2 - np.mean(y2)) ** 2)

r = 1 - (SS_res / SS_total)
r2 = 1 - (SS2_res / SS2_total)
print("R^2 NW time complexity = ", r)
print("R^2 HB time compexity  = ", r2)

# plot memory analysis
# x-axis = average length of the two aligned sequences
# y-axis = average time of alignment (alignment performed 5 times, average time shown)
nw_mem = [(19.0, 0.0234375), (479.0, 8.0703125), (1001.5, 40.83984375), (1763.5, 118.22265625), (2229.0, 213.4375), (3000, 417.546875), (4000, 722.6328125)]
hb_mem = [(19.0, 0.0039062500), (479.0, 0.1796875000), (1001.5, 0.4648437500), (1763.5, 0.9843750000), (2229.0, 1.1875000000), (3000, 1.4101562500), (4000, 1.5078125000), (5000, 1.71875), (6000, 2.0742187500), (6403.5, 2.2031250000)]
x_mem, y_mem = zip(*nw_mem)
x2_mem, y2_mem = zip(*hb_mem)

# fit quadratic polynomial to needleman
coeffs_nw_mem = np.polyfit(x_mem, y_mem, 2)
poly_nw_mem = np.poly1d(coeffs_nw_mem)
print(coeffs_nw_mem)
print(f'NW memory quadratic fit coefficients: {coeffs_nw_mem}')

# fit linear regression to hirschberg
coeffs_hb_mem = np.polyfit(x2_mem, y2_mem, 1)
poly_hb_mem = np.poly1d(coeffs_hb_mem)
print(coeffs_hb_mem)
print(f'Hirschberg memory linear fit coefficients: {coeffs_hb_mem}')

# x-vals for polynomial plotting
x_fit_mem = np.linspace(min(min(x_mem), min(x2_mem)), max(max(x_mem), max(x2_mem)), 100)
x_fit_mem_NW_only = np.linspace(min(x_mem),max(x_mem), 100)

# plotting

# needleman only
plt.figure(2)
plt.scatter(x_mem, y_mem, color='blue', label='Needleman-Wunsch')
plt.plot(x_fit_mem_NW_only, poly_nw_mem(x_fit_mem_NW_only), color='blue', linestyle='--', label = 'Quadratic Fit for Needleman-Wunsch')
plt.xlabel('Average Sequence Length')
plt.ylabel('Average Alignment Memory (MiB)')
plt.legend()

# hirschberg only
plt.figure(3)
plt.scatter(x2_mem, y2_mem, color='red', label='Hirschberg')
plt.plot(x_fit_mem, poly_hb_mem(x_fit_mem), color='red', linestyle='--', label = 'Linear Fit for Hirschberg')
plt.xlabel('Average Sequence Length')
plt.ylabel('Average Alignment Memory (MiB)')
plt.legend()

# calculating R^2 values
y_pred_mem = poly_nw_mem(x_mem)
y2_pred_mem = poly_hb_mem(x2_mem)
SS_res_mem = np.sum((y_mem - y_pred_mem) ** 2)
SS2_res_mem = np.sum((y2_mem - y2_pred_mem) ** 2)
SS_total_mem = np.sum((y_mem - np.mean(y_mem)) ** 2)
SS2_total_mem = np.sum((y2_mem - np.mean(y2_mem)) ** 2)

r_mem = 1 - (SS_res_mem / SS_total_mem)
r2_mem = 1 - (SS2_res_mem / SS2_total_mem)
print("R^2 NW memory = ", r_mem)
print("R^2 HB memory = ", r2_mem)

plt.show()