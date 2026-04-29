% Replace your orbital element extraction loop (inside the case loop) with this.
% It catches and reports any rv2coe failures, and skips bad pieces silently
% in subsequent runs.

a_hist     = nan(N, M);
e_hist     = nan(N, M);
Omega_hist = nan(N, M);
zp_hist    = nan(N, M);
za_hist    = nan(N, M);
tau_hist   = nan(N, M);

n_failed = 0;
first_failure_reported = false;

for k = 1:N
    for m = 1:M
        rk = r_traj(:,k,m);
        vk = v_traj(:,k,m);
        if any(isnan(rk)) || any(isnan(vk));  continue;  end

        try
            [aa, ee, OO, ~, ~, ~] = rv2coe(rk, vk);

            % Verify scalars
            if ~isscalar(aa) || ~isscalar(ee) || ~isscalar(OO)
                if ~first_failure_reported
                    fprintf('  [WARN] rv2coe returned non-scalar for piece %d, epoch %d\n', k, m);
                    fprintf('         size(a)=%s, size(e)=%s, size(Om)=%s\n', ...
                        mat2str(size(aa)), mat2str(size(ee)), mat2str(size(OO)));
                    first_failure_reported = true;
                end
                n_failed = n_failed + 1;
                continue
            end

            if ~isfinite(aa) || aa <= 0 || ee >= 1;  continue;  end

            a_hist(k,m)     = aa;
            e_hist(k,m)     = ee;
            Omega_hist(k,m) = OO;
            zp_hist(k,m)    = aa*(1-ee) - rE;
            za_hist(k,m)    = aa*(1+ee) - rE;
            tau_hist(k,m)   = 2*pi*sqrt(aa^3/muE) / 3600;

        catch ME
            if ~first_failure_reported
                fprintf('  [WARN] rv2coe threw at piece %d, epoch %d: %s\n', k, m, ME.message);
                first_failure_reported = true;
            end
            n_failed = n_failed + 1;
        end
    end
end
if n_failed > 0
    fprintf('  Skipped %d bad rv2coe results out of %d total points\n', n_failed, N*M);
end