# Traf 관련 파일 정리

### `additional_wireless_links_rev.m`
- MATLAB 코드 수정본  
- Adjacency Matrix(Topology)만 저장하도록 구성

***

### `graph_topologies_only.mat`
- MATLAB에서 생성한 Adjacency Matrix 저장 파일

***

### `graph_topology_test.ipynb`
- 저장된 adjacency matrix 확인용 Jupyter Notebook

***

### `trafy_test.ipynb`
- TrafPy 기능 테스트 Notebook  
- 참고 문서: https://trafpy.readthedocs.io/en/latest/tutorial_generator.html

***

### `flows_example.ipynb`
- Traf 기반 flow generation 실험 Notebook  
- 설정 (Cell 3 참고):
  - Uniform node selection  
  - 3 종류의 flow size  
  - Fixed interarrival time
 
### `flows_output_64_uniform.csv`
- 위의 flows data csv 파일  

***

### `DCN3_traf.ipynb`
- 논문 파라미터 + 현재 topology 사용하여 flow generation 수행

### `flows_output_DCN3.csv`
- 위의 flows data csv 파일
