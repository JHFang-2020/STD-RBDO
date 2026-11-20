function w = c36_wamp(S,gamma)
%C36_WAMP  Amplitude function of the non-stationary discrepancy kernel.
w = sqrt(1 + gamma*(S(:,1).^2 + S(:,2).^2));
end
