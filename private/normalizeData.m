function ndat = normalizeData(data, nAttributes)

% Author(s): Paul Blomstedt


ndat = data;
for i=1:nAttributes
    ind = find(ndat(:,i)~=0);
    if length(ind)==1
        ndat(ind,i) = 0.0001; % distinguish these from "true" zeros
    elseif length(ind)>1
        m = mean(ndat(ind,i));
        v = var(ndat(ind,i));
        if v==0 % avoid dividing by 0
            v = 1;
        end
        ndat(ind,i) = (ndat(ind,i)-m)./sqrt(v);
    end
end
