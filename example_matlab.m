% example
% Example of the bidiagonal updating algorithms: bgu and bhu 
%
% This script computes the bidiagonal factorization of B+w*p'
% The result is a new factorization Q*B1*P' in which the orthogonal
% Q and P' are stored in factored form. Depending on the algorithm these
% dense matrices may be formed explicitly.
%
%--------------------------------------------------------------------------
% 04/08/25, J.B., initial version
% 08/26/25, J.B., preparation for release
% 08/27/25, J.B., further prep. for release

clc;
clear;

% load interface
loadlibs;

% Setup problem
fc      = 200;  % 1
m       = fc*10;
n       = fc*4;
B       = [diag(1:n) + diag(2:n,1); zeros(m-n,n)];
w       = zeros(m,1);  
p       = zeros(n,1);  
w(2)    = 1;
w(6)    = 10;
p(2)    = 2;
p(3)    = 3;

if m >= 500
    w(500)  = 1;
end

rr = matlabRelease;

fprintf('***************** Bidiag update ******************* \n');
fprintf('*         Algorithms: bgu                           \n');
fprintf('*                     bhu                           \n');
fprintf('*                                                   \n');
fprintf('*           Software: Matlab %s                     \n',rr.Release);
fprintf('*             System: %s                            \n',computer);
fprintf('*                                                   \n');
fprintf('*          Prob size: m=%i                          \n',m);
fprintf('*                     n=%i                          \n',n);
fprintf('*                     nnz(w)=%i                     \n',nnz(w));
fprintf('*                     nnz(p)=%i                     \n',nnz(p));
fprintf('*                                                   \n');
fprintf('*         Release Aug. 2025                         \n');
fprintf('*         J.J. Brust (johannesbrust@yahoo.com)      \n');
fprintf('*************************************************** \n');
fprintf('\n');

% call bgu algorithm (Given's based)
tbgu = tic;
[B1,Q1,Q2,Q3,P1,P2] = bidiagup(B,w,p);
tbgu = toc(tbgu);

%%Alternatively, compute orthogonal matrices
%Im      = eye(m);
%In      = eye(n);
%%Permutation infos stores in vectors pq and pp.
% The algorithm computes Q(pq,:)*(B+w*p')*P(:,pp) = B1
%[QT,pq] = bdup_mulq(Im,Q1,Q2,Q3);
%[PT,pp] = bdup_mulp(In,P1,P2);
%Q       = QT(pq,:)';
%P       = PT(:,pp)';
%%factorization error
%err = norm(Q*B1*P-(B+w*p'),'fro');

% Reconstruct the original matrix or compare norms only (if problems are
% large)
if m <= 1000
    tmbgu   = tic;
    [BU,pq] = bdup_mulq(B1,Q1,Q2,Q3,1);
    [BU,pp] = bdup_mulp(BU,P1,P2,1);
    tmbgu   = toc(tmbgu);
    err     = norm(BU(pq,pp)-(B+w*p'),'fro');
else
    tmbgu   = tbgu; 
    err     = abs(norm(B1,'fro')-norm((B+w*p'),'fro'));
end

% Call bhu algorithm (Householder based)
% The bidiagonal elements are in the columns of B1 and Y, W
% store the essential information of Householder reflectors
% This method first needs to extract the columns, too
tbhu        = tic;
BUP         = [diag(B(1:n,1:n)),[diag(B(1:n,1:n),1);0]];
[B1h,Y,W]   = bhu(BUP,w,p,n);
tbhu        = toc(tbhu);

% Explicitly form the orthogonal matrices using the compact representation
if m <= 1000
    tmbhu   = tic;
    YL      = tril(Y);
    WL      = tril(W);
    cls     = size(Y,2);
    T       = triu(Y(1:cls,1:cls),1) + eye(cls);
    cls     = size(W,2);
    R       = triu(W(1:cls,1:cls),1) + eye(cls);
    QQ      = eye(m) - 2 * YL *( T\YL' );
    PP      = eye(n) - 2 * WL *( R\WL' );
    tmbhu   = toc(tmbhu);
    errh    = norm(B+w*p'-QQ(:,1:n)*(diag(B1h(1:n,1))+diag(B1h(1:n-1,2),1))*PP(:,1:n)','fro');
else
    tmbhu = tbhu;
    errh  = abs(norm(B1h,'fro')-norm((B+w*p'),'fro'));
end

fprintf(' alg: bgu                         \n');
fprintf(' time update: %.4f                \n',tbgu);
fprintf(' time mult:   %.4f                \n',tmbgu);
fprintf(' error:       %.5e                \n',err);
fprintf('                                  \n');
fprintf(' alg: bhu                         \n');
fprintf(' time update: %.4f                \n',tbhu);
fprintf(' time mult:   %.4f                \n',tmbhu);
fprintf(' error:       %.5e                \n',errh);

