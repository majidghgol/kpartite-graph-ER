function [B, C, fs ] = WeightRefinement(B, C, A, r, powerA, res, hasCon, w)
%  
%% iteration
numchar=0; 
fs=zeros(1,res.numit);
myzero= 1e-9;
k=size(A,1);
for it=1:res.numit
    for i=1:k
        N=zeros(size(C{i}));
        D=zeros(size(C{i}));
        for j=i:k
            if hasCon(i,j)
                if j~=i
                    N=N+w(i,j)*A{i,j}*C{j}*B{i,j}';
                    D=D+w(i,j)*C{i}*B{i,j}*C{j}'*C{j}*B{i,j}';
                else
                    N=N+w(i,j)*A{i,j}*C{i};
                    D=D+w(i,j)*diag(sum(A{i,j}))*C{i};
                end
            end
        end
        C{i}=C{i}.*sqrt(N./max(D,myzero));
        %s=cputime;
        C{i}=C{i}./repmat(sum(C{i},2)+myzero,1,r(i));
%         t=cputime;
%         disp('normalization time');
%         disp(t-s);
    end
    
    for i=1:k
        for j=i+1:k % save speed do half the updates, implies SYMMETRIC graph i.e. nondirected
            if hasCon(i,j)
				if j~=i
					B{i,j}=B{i,j}.*sqrt((C{i}'*A{i,j}*C{j})./max(C{i}'*C{i}*B{i,j}*C{j}'*C{j},myzero));
					B{j,i}=B{i,j}';
				end
            end
        end
    end
    f=0;
    for i=1:k
        for j=i:k % save speed do half the updates, implies SYMMETRIC graph i.e. nondirected
            if hasCon(i,j)
                if j~=i
                    f=f+w(i,j)*sum(sum((A{i,j}-C{i}*B{i,j}*C{j}').^2));
                else
                    f=f+w(i,j)*sum(sum((A{i,i}-C{i}*C{i}').^2));
                end
            end
        end
    end
    fs(it)=f;
    if res.verbose
        if res.verbosecompact 
            for i=1:numchar; fprintf('\b'); end
        else
            fprintf('\n');
        end
    	s=sprintf('it=%i  cost=%8.4f',it,f);
        numchar=length(s);fprintf(s);
    end
    
     if it>1 && (fs(it-1)-fs(it))/powerA<res.abort
         if res.verbose, fprintf('\ncovergence after %i iterations\n',it); end
         break;
     end
end

