ccc

%create problem:
N = 10;
lb = repmat([-0.5;0],N/2,1);
ub = repmat([inf;0.5],N/2,1);
obj = qpBoxFcnTest('N',N,...
    'lb',lb,...
    'ub',ub,...
    'tolFun',1e-8,'maxIter',50);
obj.setup(1,1)

opts = optimoptions(@quadprog,'Display','None',...
    'Algorithm','trust-region-reflective');

iter = 50;
T = zeros(iter,1);
err = T;
errCount = 0;
victCount = 0;
for n = 1:iter
    a = randi([2 10]);
    b = randn;
    tic
    [x1,chk] = obj.step(a,b);
    if chk ~= 2
        warning('ran out of iterations')
    end
    dt = toc;
    H = obj.hCalc(a,b);    
    f = obj.fCalc(a,b);
    q1 = obj.fQuad(H,f,x1);
    
    %compare against quadprog
    tic
    [x2,q2] = quadprog(2*H,f,[],[],[],[],lb,ub,[],opts);
    dt2 = toc;
    err(n) = 100*(q1-q2)/q2;
    T(n) = dt2/dt;
    fprintf('\n%0.2f times faster\nError: %0.5f%%\n',T(n),err(n))
    if q1 < q2
        victCount = victCount + 1;
    end
    if err(n) < -1 && q1 > obj.tolFun
        errCount = errCount + 1;
        disp([q1 q2])
        warning('qpBox didn''t perform too well here')
    end
end
fprintf('\nNumber of poor performances: %d\n',errCount)
fprintf('\nNumber of times qpBox beat quadprog: %d\n',victCount)
fprintf('Average error: %0.2e%%\n',abs(mean(err(err < 0))))
fprintf('Average execution time ratio: %0.2f\n',mean(T))