function flag = triadCondition(q,p,k)

flag = 0;
if (q >= (k-p) )
    if ( q <= (k+p) )
        if (q >= (p-k) )
            flag = 1;
        end
    end
end

end