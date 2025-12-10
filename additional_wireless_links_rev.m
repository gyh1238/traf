%%===============================
%%  Wireless Torus: Topology Only
%%  - G_ori  : torus graph (unweighted, undirected)
%%  - G4_add : torus + additional links (unweighted, undirected)
%%  - G5_add : same topology with degree-based weights (digraph)
%%===============================

clear; clc;

%% ---- Parameters ----
n_size_max = 8;
n_size_min = 4;
n_loop     = 2;                    % how many graphs per (size, del_per, pa)

del_per_list = 10:20:150;          % percentage of additional links
per_size     = numel(del_per_list);

pa_list = 0:1/3:1;                 % preferential attachment ratio
n_pa    = numel(pa_list);

% (원코드에는 pa_cut_th1, pa_cut_th2 있었지만 결국 nn_size^3으로 덮어씀)

%% ---- Storage for adjacency matrices ----
graph_k4_add   = zeros(4^3,4^3,n_loop,per_size,n_pa);
graph_k6_add   = zeros(6^3,6^3,n_loop,per_size,n_pa);
graph_k8_add   = zeros(8^3,8^3,n_loop,per_size,n_pa);

graph_k4_torus = zeros(4^3,4^3,n_loop,per_size,n_pa);
graph_k6_torus = zeros(6^3,6^3,n_loop,per_size,n_pa);
graph_k8_torus = zeros(8^3,8^3,n_loop,per_size,n_pa);

%% ---- Main loop over torus size ----
for nn_size = n_size_min:2:n_size_max      % 4, 6, 8
    
    % pa 최대 degree 비율을 의미하던 변수지만,
    % 원코드에서 바로 아래에서 nn_size^3으로 덮어씀
    pa_cut_th1 = 0.1;
    pa_cut_th2 = 0.03;
    if nn_size < 6
        pa_cut_num = ceil(nn_size^3 * pa_cut_th1);
    else
        pa_cut_num = ceil(nn_size^3 * pa_cut_th2);
    end
    % 실제 사용되는 건 이 값 (degree 제한 사실상 해제)
    pa_cut_num = nn_size^3;
    
    %% ---- 1) 3D torus 좌표 및 기본 adjacency (graph_ori, G_ori) ----
    hpc_size      = nn_size;
    hpc_cor_temp  = 0:hpc_size^3-1;
    hpc_cor_temp2 = dec2base(hpc_cor_temp, hpc_size);
    hpc_cor       = zeros(size(hpc_cor_temp2,1), 3);
    
    for ii = 1:size(hpc_cor_temp2,1)
        val1 = base2dec(hpc_cor_temp2(ii,1), hpc_size);
        val2 = base2dec(hpc_cor_temp2(ii,2), hpc_size);
        val3 = base2dec(hpc_cor_temp2(ii,3), hpc_size);
        hpc_cor(ii,:) = [val1, val2, val3];
    end
    
    % Euclidean / Hamming distance
    dist_eu  = pdist2(hpc_cor, hpc_cor);
    dist_ham = pdist2(hpc_cor, hpc_cor, 'hamming');
    
    % 인접 노드 (거리 1) edge
    dist_edge = double(abs(dist_eu - 1) < 0.01);
    
    % Torus wrap-around edge
    zz  = abs(dist_eu  - (hpc_size - 1)) < 0.001;
    zz2 = abs(dist_ham - 1/3)            < 0.001;
    torus_add        = zz .* zz2;
    dist_edge_torus  = dist_edge + torus_add;
    
    graph_ori = dist_edge_torus;
    G_ori     = graph(graph_ori);         % unweighted torus (undirected)
    
    %% ---- 2) Wireless conversion (추가 링크) ----
    for nn_per = 1:per_size               % 추가 링크 비율 루프
        del_per = del_per_list(nn_per);   % [%]
        
        for nn_pa = 1:n_pa                % preferential attachment 비율 루프
            
            % nn_loop 개의 valid 그래프를 뽑을 때까지 반복
            loop = 1;
            while loop <= n_loop
                
                % 기준: torus adjacency에서 시작
                graph_pa_temp2 = dist_edge_torus;
                node_pa_temp2  = size(graph_pa_temp2,1);
                
                % 현재 torus에서 존재하는 edge 수
                low_tri    = tril(dist_edge_torus,-1);
                ind_tri    = low_tri(:);
                con_loc_before = find(abs(ind_tri - 1) < 0.001);
                
                % 추가할 edge 개수 (퍼센트 기반)
                del_num = ceil(del_per * numel(con_loc_before) / 100);
                
                cutoff_th = 0.7;  % early 단계에는 degree^2 weight 사용
                
                %% --- 추가 링크 붙이기 (preferential / random) ---
                for pa_num = 1:del_num
                    G_temp2 = graph(graph_pa_temp2);
                    
                    if rand(1) < pa_list(nn_pa)
                        % ----- preferential attachment -----
                        r_check2 = 1;
                        while r_check2 > 0
                            d_temp2    = degree(G_temp2);
                            d_min2     = min(d_temp2);
                            d_min_loc2 = find(abs(d_temp2 - d_min2) < 0.001);
                            
                            % degree 가장 낮은 노드 중 하나 선택
                            cor_i2 = datasample(d_min_loc2, 1);
                            
                            % 나머지 한 노드는 degree-based weight로 선택
                            if pa_num < ceil(del_num * cutoff_th)
                                cor_j2 = datasample(1:node_pa_temp2, 1, ...
                                                    'Weights', degree(G_temp2).^2);
                            else
                                cor_j2 = datasample(1:node_pa_temp2, 1, ...
                                                    'Weights', degree(G_temp2));
                            end
                            
                            degree_chk = degree(G_temp2);
                            if graph_pa_temp2(cor_i2,cor_j2) == 0 && ...
                               cor_i2 ~= cor_j2 && ...
                               degree_chk(cor_j2) < pa_cut_num
                                
                                graph_pa_temp2(cor_i2,cor_j2) = 1;
                                graph_pa_temp2(cor_j2,cor_i2) = 1;
                                r_check2 = -10;
                            end
                        end
                        
                    else
                        % ----- random attachment -----
                        rr_check2 = 1;
                        while rr_check2 > 0
                            cor_i2 = datasample(1:node_pa_temp2,1);
                            cor_j2 = datasample(1:node_pa_temp2,1);
                            
                            if graph_pa_temp2(cor_i2,cor_j2) == 0 && ...
                               cor_i2 ~= cor_j2
                                graph_pa_temp2(cor_i2,cor_j2) = 1;
                                graph_pa_temp2(cor_j2,cor_i2) = 1;
                                rr_check2 = -10;
                            end
                        end
                    end
                end % end for pa_num
                
                % 최종 추가링크 포함 adjacency
                graph_pa_add = graph_pa_temp2;
                G4_add       = graph(graph_pa_add);   % unweighted, undirected
                
                %% --- 3) degree 기반 weight 붙인 digraph (G5_add) ---
                graph_pa_weight_add = (~~graph_pa_add .* degree(G4_add)).';
                G5_add              = digraph(graph_pa_add);
                G5_add.Edges.Weight = nonzeros(graph_pa_weight_add);
                
                % distance matrix → 연결성 체크용 (원래 latency 계산용이었음)
                distance_add = distances(G5_add);
                
                % 그래프가 완전히 끊어지지 않았을 때만 저장
                if ~isinf(mean(distance_add(:)))
                    
                    if nn_size == 4
                        graph_k4_add(:,:,loop,nn_per,nn_pa)   = graph_pa_add;
                        graph_k4_torus(:,:,loop,nn_per,nn_pa) = graph_ori;
                    elseif nn_size == 6
                        graph_k6_add(:,:,loop,nn_per,nn_pa)   = graph_pa_add;
                        graph_k6_torus(:,:,loop,nn_per,nn_pa) = graph_ori;
                    else % nn_size == 8
                        graph_k8_add(:,:,loop,nn_per,nn_pa)   = graph_pa_add;
                        graph_k8_torus(:,:,loop,nn_per,nn_pa) = graph_ori;
                    end
                    
                    if mod(loop,10) == 1   % 디버깅용 출력
                        [nn_size, nn_per, nn_pa, loop]
                    end
                    
                    loop = loop + 1;
                end
                
            end % while loop <= n_loop
            
        end % nn_pa
    end % nn_per
    
end % nn_size

%% 필요하면 여기서 저장
save('graph_topologies_only.mat', ...
     'graph_k4_add','graph_k6_add','graph_k8_add', ...
     'graph_k4_torus','graph_k6_torus','graph_k8_torus');
