function [optX,optFcn,noIter] = GSO(objf,solSize,userOptions)
% Programmed by Daniel Alienyi

if nargin< 2
    optFcn = [];    %empty matrix optFcn, ie it is a zero element, no result 
    optX = [];      %empty matrix optX, ie it is a zero element, no result 
    noIter = [];     % empty matrix noIter, ie it is a zero element, no result
    return     %terminate
else
    % Initialising the options
    options.population = 100;
    options.L0 = 5;
    options.r0 = 3;
    options.rho = 0.4;
    options.y = 0.6;
    options.B = 0.08;
    options.s = 0.6;
    options.rs = 10;
    options.nt = 10;
    options.maxIter = 1000;
    options.xmin = -20;
    options.xmax = 20;
    
    options.showplot = true;
    
    if nargin > 2
        allFields = fieldnames(userOptions);
        fdSz = length(allFields);
        for j = 1:fdSz;
            options.(allFields{j}) = userOptions.(allFields{j});
        end
    end
end

noPop = options.population;
L0 = options.L0;
r0 = options.r0;
rho = options.rho;
y = options.y;
B = options.B;
s = options.s;
rs = options.rs;
nt = options.nt;
maxIter = options.maxIter;
xmin = options.xmin;
xmax = options.xmax;

objfcn = @(x)ConvertToMin(x, objf);

Xs = (xmax - xmin)*rand(noPop,solSize) + xmin;  % Initialise Xs. its xi in paper

luciferin = L0*ones(noPop,1);    % Initialising the luciferin
decision_range = r0*ones(noPop,1);  % Initialising the decision range
numList = 1:noPop;
iter = 1;

%% ---------- Start of plot --------------------------------------------
if solSize < 3 & options.showplot     % Determines if a plot should be made
    fig = figure;
end
if solSize == 1 & options.showplot
    oneDimPlotStarter(objfcn,xmin,xmax)
elseif solSize == 2
    twoDimPlotStarter(objfcn,xmin,xmax)
end
if solSize < 3 & options.showplot
    pause(0.2)
    hold on
end
if solSize == 1
    hd = oneDimPlotUpdate(objfcn,Xs);
elseif solSize == 2
    hd = twoDimPlotUpdate(objfcn,Xs);
end
%% ------- End of plot -----------------------------------------------------
while iter <= maxIter
    
    if solSize < 3 & options.showplot % Updates the plot
        pause(0.2);
        delete(hd)
    end
    
    % Updating the luciferin
    luciferin = (1-rho)*luciferin + y*getFcn(objfcn,Xs);
    
    [bestL,bestPos] = max(luciferin);
    
    % Moving the Glow-worms
    for ii = 1:noPop;
        curX = Xs(ii,:);
        curLuciferin = luciferin(ii);
        distFromI = EuclidDistance(Xs,repmat(curX,noPop,1));
        Ni = find((distFromI < decision_range(ii)) & (luciferin > curLuciferin));
        if isempty(Ni)  % If no glow-worm exists within its local range
            Xs(ii,:) = curX;
        else
            localRangeL = luciferin(Ni);
            localRangeX = Xs(Ni,:);
         
            probs = (localRangeL - curLuciferin)/sum(localRangeL - curLuciferin);
            selectedPos = SelectByRoulete(probs);
            selectedX = localRangeX(selectedPos,:);
            Xs(ii,:) = curX + s*(selectedX - curX)/EuclidDistance(selectedX,curX);
        end
        neighborSz = length(Ni);
        decision_range(ii) = min([rs,max([0,decision_range(ii) + B*(nt-neighborSz)])]);
    end
    
    iter = iter + 1;
    
    % Updating the plot
    if solSize == 1 & options.showplot
        hd = oneDimPlotUpdate(objfcn,Xs);
    elseif solSize == 2
        hd = twoDimPlotUpdate(objfcn,Xs);
    end
    
end

if solSize < 3 & options.showplot
    hold off
end

optX = Xs(bestPos,:);
optFcn = getFcn(objfcn,optX);
noIter = iter - 1;


function Fx = getFcn(objfcn,Xs)
n = size(Xs,1);
Fx = ones(n,1);
for k = 1:n;
    Fx(k) = objfcn(Xs(k,:));
end

function ret = EuclidDistance(pos1,pos2)
ret = sqrt(sum((pos1-pos2).^2,2));

function ret = SelectByRoulete(allProb)
cumProb = cumsum(allProb);
rn = rand;
hd = find(cumProb >= rn);
ret = hd(1);

function oneDimPlotStarter(objfcn,mn,mx)
    xrange = linspace(mn,mx,500);
    yrange = objfcn(xrange);
    plot(xrange,yrange);
    
function hd = oneDimPlotUpdate(objfcn,X)
    hd = plot(X,objfcn(X),'ro','markerfacecolor','m');
    
function twoDimPlotStarter(objfcn,mn,mx)
    rg = linspace(mn,mx,120);
    [xrange,yrange] = meshgrid(rg,rg);
    zrange = xrange;
    sz = length(rg);
    for x = 1:sz;
        for y = 1:sz;
            zrange(x,y) = objfcn([xrange(x,y),yrange(x,y)]);
        end
    end
    surf(xrange,yrange,zrange);
    shading(gca,'interp')
    %contour(xrange,yrange,zrange);
    
function hd = twoDimPlotUpdate(objfcn,X)
    zrange = getFcn(objfcn,X);
    xrange = X(:,1);
    yrange = X(:,2);
    hd = plot3(xrange,yrange,zrange,'ro','markersize',8,'markerfacecolor','m');
    
function hd = twoDimNeighborPlot(objfcn,X, probs)
    zrange = getFcn(objfcn,X);
    xrange = X(:,1);
    yrange = X(:,2);
    mn = min(probs);
    mz = max(probs);
    largeness = 9 + 6*(probs - mn)./(mz - mn);
    sz = length(zrange);
    hd = zeros(1,sz);
    for i = 1:sz
        hd = [hd,plot3(xrange(i),yrange(i),zrange(i),'go','markersize',largeness(i),'markerfacecolor','g')];
    end
    

function minObjFcn = ConvertToMin(x, objfcn)
    fcn = objfcn(x);
    if fcn >= 0
        minObjFcn = 1/(1+fcn);
    else
        minObjFcn = 1 + abs(fcn);
    end
  