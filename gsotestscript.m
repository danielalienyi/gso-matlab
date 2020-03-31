fcn1 = @(x)0.5 + ((sin(sqrt(sum(x.^2)))).^2 - 0.5)./(1 + 0.001*sum(x.^2)).^2;

fcn2 = @(x)sum(x.^2);

fcn3 = @(x)sum(x.^2 - 10*cos(2*pi*x) + 10);

fcn4 = @(x)1 + sum(x.^2)/4000 - prod(cos(x./(1:length(x))));

fcn5 = @(x)sum(0.2*x.^2 + 0.1*x.^2.*sin(2*x));

fcn6 = @(x)sum(x.^2)^0.25*(sin(50*sum(x.^2)^0.1)^2 + 1.0);


opt1 = struct('xmin', -20, 'xmax', 20, 'r0', 5, 'rs', 10, 'population', 500, 's', 0.5,'maxIter', 100, 'showplot', 1);
opt2 = struct('xmin', -20, 'xmax', 20, 'population', 500,'maxIter', 200, 'showplot', 1);
opt3 = struct('xmin', -7, 'xmax', 7, 'population', 500,'maxIter', 100, 'r0', 0.5, 'rs', 3, 's', 0.001, 'showplot', 1);
opt4 = struct('xmin', -20, 'xmax', 20, 'population', 500,'maxIter', 300, 'r0', 2, 'rs', 5, 'showplot', 1);
opt5 = struct('xmin', -8, 'xmax', 8, 'population', 500,'maxIter', 100, 'r0', 1, 'rs', 5, 's', 0.1, 'showplot', 1);
opt6 = struct('xmin', -20, 'xmax', 20, 'population', 500,'maxIter', 200);

GSO(fcn2, 2, opt2);
GSO(fcn4, 2, opt4);
GSO(fcn5, 2, opt5);
