"""
Tests for the class Simulation.

First version by R. Ferry on February 2021.
"""
import numpy as np
import sys
sys.path.append("../PythonFastCycles/")
from simulation import Simulation
import pytest

# def test_create_all_files(tmp_path):
#     Test = Simulation(str(tmp_path), 'Test', 30*1e9, 0.0075, 0.01)
#     sigma_dot = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.01], [0.0, 0.01, 0.0]])
#     Test.create_all_files(2, sigma_dot, geom_type='1fault')
    
#     with open(str(tmp_path) + 'Test/config.in', 'r') as file:
#         content = file.readlines()
#     assert content[0]=='&main\n'

    
# class TestCreateAllFiles():
#     @pytest.fixture(scope="class")
#     def setup_class(self, tmp_path):
#         Test = Simulation(str(tmp_path), 'Test', 30*1e9, 0.0075, 0.01)
#         sigma_dot = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.01], [0.0, 0.01, 0.0]])
#         Test.create_all_files(2, sigma_dot, geom_type='1fault') 
        
#         self.tmp_path = tmp_path

#     @pytest.fixture 
#     def test_create_config(self):
#         with open(self.tmp_path + 'Test/config.in', 'r') as file:
#             content = file.readlines()
#         assert content[0]=='&main'

@pytest.fixture(scope="class")
def setup(tmp_path):
    Test = Simulation(str(tmp_path), 'Test', 30*1e9, 0.0075, 0.01)
    sigma_dot = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.01], [0.0, 0.01, 0.0]])
    Test.create_all_files(2, sigma_dot, geom_type='1fault') 
        
    with open(str(tmp_path) + 'Test/config.in', 'r') as file:
        content = file.readlines()
    
    return content

class TestCreateAllFiles():
    def test_create_config(self, setup):
        assert content[0]=='&main'

    
