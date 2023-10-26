function evalres = brief_display_solution(asol,gt,OK);

if nargin<3,
    %OK = ones(size(sgt_resamp(1,:)));
    OK = isfinite(gt.sgt_resamp(1,:));
end

evalres = evaluate_solution(asol,gt,OK);
display(['#ok s-punkter ' num2str(sum(evalres.sok.*OK)) '/' num2str(sum(OK))]);

