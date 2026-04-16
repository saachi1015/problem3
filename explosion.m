function [v, bc] = explosion(mu_v, sigma_v, mu_bc, sigma_bc, N)
v = mu_v.' + sigma_v.* randn(N,1);
v = v.';
pd_ball = makedist('Lognormal','mu', mu_bc,'sigma', sigma_bc);
bc = random(pd_ball, N, 1);
bc = bc.';
end