f = @(x) (x-2).^(1/4);

% Level 1
display 'Level 1'
Si1 = simpson(f,3,5,1)
S11 = simpson(f,3,4,1)
S21 = simpson(f,4,5,1)
St1 = S11 + S21
err = abs(St1 - Si1)/15

% Level 2
display 'Level 2'
Si2 = S11
S12 = simpson(f,3,3.5,1)
S22 = simpson(f,3.5,4,1)
St2 = S12 + S22
err = abs(St2 - Si2)/15

% Level 2
display 'Level 2'
Si3 = simpson(f,4,5,1)
S13 = simpson(f,4,4.5,1)
S23 = simpson(f,4.5,5,1)
St3 = S13 + S23
err = abs(St3 - Si3)/15


