function power_intencity_plot()
Power_int=importdata('D:\h04202022\intensity over power.txt');
Power_int=Power_int.data;
Power_int_mean=Power_int(:,[2 4 6]);
Power_int_std=Power_int(:,[3 5 7]);
plot(Power_int_mean)
hold on
errorbar(Power_int_mean,Power_int_std)
x
