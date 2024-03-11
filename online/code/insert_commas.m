function comma_N = insert_commas(N)
% Insert commas into numbers, e.g. 1234567 -> '1,234,567'.

N = num2str(N);
num_digits = length(N);
comma_N = '';

for j = 1:num_digits
    comma_N = strcat(N(num_digits-j+1),comma_N);
    if ~mod(j,3)
        comma_N = strcat(',',comma_N);
    end    
end
if strcmp(comma_N(1), ',')
    comma_N(1) = [];
end