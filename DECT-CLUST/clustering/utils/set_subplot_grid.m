function ha = set_subplot_grid(K)
if(K==2)
    ha = tight_subplot(1,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==3)
    ha = tight_subplot(1,3, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==4)
    ha = tight_subplot(2,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==5||K==6)
    ha = tight_subplot(3,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==7||K==8)
    ha = tight_subplot(4,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==9)
    ha = tight_subplot(3,3, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>9 && K<=12)
    ha = tight_subplot(4,3, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>12 && K<=16)
    ha = tight_subplot(4,4, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>6 && K<=20)
    ha = tight_subplot(4,5, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>20 && K<=25)
    ha = tight_subplot(5,5, [0.02 0.06]);
elseif (K>25 && K<=30)
    ha = tight_subplot(6,5, [0.02 0.06]);
elseif (K>30 && K<=36)
    ha = tight_subplot(6,6, [0.02 0.06]);
elseif (K>36 && K<=42)
    ha = tight_subplot(7,6, [0.02 0.06]);
elseif (K>35 && K<=48)
    ha = tight_subplot(8,6, [0.02 0.06]);
elseif (K>48 && K<=54)
    ha = tight_subplot(9,6, [0.02 0.06]);
elseif (K>54 && K<=63)
    ha = tight_subplot(9,7, [0.02 0.06]);
else
    warning('Please add the plot case for this many clusters');
end

end