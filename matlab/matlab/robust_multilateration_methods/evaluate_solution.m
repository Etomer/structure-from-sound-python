function evalres = evaluate_solution(asol,gt,OK);

if nargin<3,
    %OK = ones(size(sgt_resamp(1,:)));

    OK = isfinite(gt.sgt_resamp(1,:));
end

threshold = 0.15;

% Check how many measurements there are for each time instant
tn = size(asol.s,2);
hh = hist(asol.J,1:tn);

% Do procrustes first
[D, Z, TRANSFORM] = procrustes(gt.rgt',asol.r','Scaling',false);
r_new = (TRANSFORM.b * asol.r' * TRANSFORM.T + TRANSFORM.c)';
s_new = (TRANSFORM.b * asol.s' * TRANSFORM.T + repmat(TRANSFORM.c(1,:),size(asol.s,2),1))';

sdiff = sqrt(sum( (gt.sgt_resamp - s_new).^2 ));
sok = sdiff<threshold;

evalres.hh = hh;
evalres.sok = sok;

tmp = sok & (hh>=3);
tmp = sok;

evalres.nrok = sum(tmp.*OK);
evalres.maxok = sum(OK);
evalres.r_new = r_new;
evalres.s_new = s_new;


%%

