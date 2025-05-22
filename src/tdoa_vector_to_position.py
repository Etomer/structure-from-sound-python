import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

import json
import subprocess




def matlab_tdoa_vector_to_positions(input_folder, output_folder=None):
    config = json.load(open("config.json","r"))
    if output_folder == None:
        output_folder = input_folder

    # resaves the tdoa_vectors to csv so matlab can read them
    tdoav = np.load(os.path.join(input_folder,"tdoa_vectors.npy"))
    df = pd.DataFrame(tdoav)
    df.to_csv(os.path.join(output_folder ,"tdoa_vectors_to_matlab.csv"))

    # calls the matlab code
    a = subprocess.run(config["matlab"] + " -nodesktop -nosplash -nodisplay -r " + "\" addpath('./matlab/matlab');tdoa_for_python('" +
            output_folder + "');exit\"", shell = True)


def python_tdoa_vector_to_positions(input_folder, output_folder=None, measurement_noise_std = 0.1):
    if output_folder == None:
        output_folder = input_folder
    
    # resaves the tdoa_vectors to csv so matlab can read them
    tdoa = np.load(os.path.join(input_folder,"tdoa_vectors.npy")).T
    tp = TxoaProblem(tdoa, tol=measurement_noise_std)
    tp.solve_for_offset(OffsetSolver95)
    #tp.bundle(steps=10000, lr=1e-3)
    tp.ransac_expand_to_all_rows()
    tp.bundle(steps=10000, lr=1e-3)
    R,S = tp.upgrade(ransac_iter=10000)

    pd.DataFrame(R).to_csv(os.path.join(output_folder ,"receiver_positions.csv"), header=False,index=False)
    pd.DataFrame(S).to_csv(os.path.join(output_folder ,"sender_positions.csv"), header=False,index=False)
    pd.DataFrame(tp.sol.o).to_csv(os.path.join(output_folder ,"offsets.csv"), header=False,index=False)


# Helper functions ------------------------------------------------------------------

class UvaboSolution():

    def __init__(self, n_receivers, n_senders, dim=3):
        u = np.empty((3, n_receivers))
        u[:] = np.nan
        self.u = u

        v = np.empty((3, n_senders))
        v[:] = np.nan
        self.v = v

        a = np.empty((n_receivers,1))
        a[:] = np.nan
        self.a = a

        b = np.empty((1, n_senders))
        b[:] = np.nan
        self.b = b

        o = np.empty((1,n_senders))
        o[:] = np.nan
        self.o = o

from enum import Enum
from abc import ABC, abstractmethod

def init_uvab(d):
    dsquare = d**2
    doubleCompaction = dsquare - dsquare[0:1,:] - dsquare[:,0:1] + dsquare[0,0]
    uu,ss,vv = np.linalg.svd(doubleCompaction/-2)
    u =uu[:,:3].T
    v = np.diag(ss[:3])@vv[:3,:]
    a = dsquare[:,0:1] - dsquare[0,0]/2
    b = dsquare[0:1,:] - dsquare[0,0]/2
    return u,v,a,b

class TxoaType(Enum):
    TOA = 1
    TDOA = 2
    COTDOA = 3

class TxoaProblem():
    def __init__(self, data, problem_type=TxoaType.TDOA, dim=3, tol = 0.1, tol_std=2.0):
        self.data = data
        self.sol = UvaboSolution(*data.shape, dim=dim)
        self.problem_type = problem_type
        self.dim = dim
        self.tol = tol
        self.tol_std = tol_std # number of standard deviations to include as inliers

    def solve_for_offset(self, solver, outer_ransac_iters=100, inner_ransac_iters=10):   
        #sample_subset = lambda : self.data[np.random.permutation(self.data.shape[0])[:solver.get_needed_receivers()],np.random.permutation(self.data.shape[1])[:solver.get_needed_senders()]]

        most_inliers = -1
        best_sol = None
        for _ in range(outer_ransac_iters):
            mics_choice = np.random.permutation(self.data.shape[0])[:solver.get_needed_receivers()]
            sound_choice = np.random.permutation(self.data.shape[1])[:solver.get_needed_senders()]
            data_subset = self.data[mics_choice][:,sound_choice]
            try:
                offsets = solver.solve(data_subset)
            except:
                continue
            d = data_subset - offsets
            u,v,a,b = init_uvab(d)

            cur_solution = UvaboSolution(*self.data.shape,dim=self.dim)

            cur_solution.u[:,mics_choice] = u
            cur_solution.v[:,sound_choice] = v
            cur_solution.a[mics_choice] = a
            cur_solution.b[:,sound_choice] = b
            cur_solution.o[:,sound_choice] = offsets
            cur_problem = TxoaProblem(self.data, problem_type=self.problem_type, dim=self.dim)
            cur_problem.sol = cur_solution

            cur_problem.ransac_expand_to_all_cols(ransac_iter=inner_ransac_iters)

            res = -2*cur_solution.u.T@cur_solution.v + cur_solution.a + cur_solution.b - (self.data - cur_solution.o)**2

            # Since the noise is in self.data, we have to rescale errors in res.
            res_normer = np.abs(2*(self.data - cur_solution.o)*self.tol) + self.tol**2
            

            if np.sum(np.abs(res/res_normer) < self.tol_std) > most_inliers:
                most_inliers = np.sum(np.abs(res/res_normer) < self.tol_std)
                best_sol = cur_solution
                compaction_row = mics_choice[0]
                compaction_col = sound_choice[0]
        print(most_inliers)
        #print(np.median(best_sol.o[0,sound_choice] - offsets_gt[0,sound_choice]))
        self.sol = best_sol
        self.compaction_row = compaction_row
        self.compaction_col = compaction_col

    def ransac_expand_col(self, new_sound_idx, ransac_iter=10):
        needed_eqs = self.dim + 2
        most_inliers = -1
        best_sol = None
        known_mics = np.argwhere(np.logical_not(np.isnan(self.sol.a)))[:,0]

        for _ in range(ransac_iter):
            mic_local_choice = np.random.permutation(known_mics.shape[0])[:needed_eqs]
            M = np.concatenate([2*self.data[known_mics[mic_local_choice],new_sound_idx:new_sound_idx+1],
                                -np.ones((needed_eqs,1)),
                                -2*self.sol.u[:,known_mics[mic_local_choice]].T,
                                np.ones((needed_eqs,1))
                                ],axis=1) # variable order is [o_j, o_j^2, v_j, b_j]
            B = self.data[known_mics[mic_local_choice],new_sound_idx:new_sound_idx+1]**2 - self.sol.a[known_mics[mic_local_choice]]
            if np.any(np.isnan(M)) or np.any(np.isnan(B)):
                print("O no! found none")
                continue
            x_part,_,_,_ = np.linalg.lstsq(M,B,rcond=None)
            x_hom = np.array([0,1,*([0]*self.dim),1])
            x_hom = np.expand_dims(x_hom,1)
            res = x_part[0]**2 - x_part[1]
            x = x_hom*res + x_part 

            o_new = x[0]
            v_new = x[2:5]
            b_new = x[5]

            lh = np.expand_dims((self.data[known_mics,new_sound_idx] - o_new)**2,axis=1)
            rh = -2*self.sol.u[:,known_mics].T@v_new+self.sol.a[known_mics]+b_new

            res = lh - rh
            # Since the noise is in self.data, we have to rescale errors in res.
            res_normer = np.abs(2*(self.data[known_mics,new_sound_idx] - o_new)*self.tol) + self.tol**2
            if np.sum(np.abs(res)/res_normer < self.tol_std) > most_inliers:
                most_inliers = np.sum(np.abs(res)/res_normer < self.tol_std)
                best_sol = (o_new[0], v_new[:,0], b_new[0])
        self.sol.o[0,new_sound_idx] = best_sol[0]
        self.sol.v[:,new_sound_idx] = best_sol[1]
        self.sol.b[0,new_sound_idx] = best_sol[2]

    def ransac_expand_row(self, new_mic_idx, ransac_iter=10):
        
        known_sounds = np.argwhere(np.logical_not(np.isnan(self.sol.b)))[:,1]
        needed_eqs = self.dim + 1
        most_inliers = -1
        best_sol = None
        for _ in range(ransac_iter):
            sound_choice = known_sounds[np.random.permutation(known_sounds.shape[0])[:needed_eqs]]
            M = np.concatenate([-2*self.sol.v[:,sound_choice].T,
                                np.ones((needed_eqs,1))
                                ],axis=1) # variable order is [o_j, o_j^2, v_j, b_j]
            B = (self.data[new_mic_idx,sound_choice] - self.sol.o[0,sound_choice])**2 - self.sol.b[0,sound_choice]
            #B = tdoa[mics_choice[mic_local_choice],new_sound_idx:new_sound_idx+1]**2 - a[mic_local_choice] 
            if np.any(np.isnan(M)) or np.any(np.isnan(B)):
                print("O no! found none")
                continue
            x = np.linalg.solve(M,B)

            u_new = np.expand_dims(x[:3],axis=1)
            a_new = x[3]

            lh = np.expand_dims((self.data[new_mic_idx,known_sounds] - self.sol.o[:,known_sounds])**2,axis=1)
            rh = -2*u_new.T@self.sol.v[:,known_sounds]+a_new+self.sol.b[:,known_sounds]

            res = lh - rh
            res_normer = np.abs(2*(self.data[new_mic_idx,known_sounds] - self.sol.o[:,known_sounds])*self.tol) + self.tol**2
            if np.sum(np.abs(res)/res_normer < self.tol_std) > most_inliers:
                most_inliers = np.sum(np.abs(res)/res_normer < self.tol_std)
                best_sol = (u_new, a_new)
        self.sol.u[:,new_mic_idx] = best_sol[0][:,0]
        self.sol.a[new_mic_idx] = best_sol[1]

    def ransac_expand_to_all_cols(self, ransac_iter=10):
        known_sounds = np.argwhere(np.logical_not(np.isnan(self.sol.b)))[:,1]
        for new_sound in np.setdiff1d(np.arange(self.data.shape[1]),known_sounds):
            self.ransac_expand_col(new_sound,ransac_iter=ransac_iter)
    
    def ransac_expand_to_all_rows(self, ransac_iter=10):
        known_mics = np.argwhere(np.logical_not(np.isnan(self.sol.a)))[:,0]
        for new_mic in np.setdiff1d(np.arange(self.data.shape[0]),known_mics):
            self.ransac_expand_row(new_mic,ransac_iter=ransac_iter)
        
    def get_residuals(self):
        return (-2*self.sol.u.T@self.sol.v + self.sol.a + self.sol.b) - (self.data - self.sol.o)**2

    def get_normed_residuals(self):
        temp = (-2*self.sol.u.T@self.sol.v + self.sol.a + self.sol.b) - (self.data - self.sol.o)**2
        normer = 2*np.abs(self.data - self.sol.o) + self.tol**2
        return temp/normer


    def bundle(self, lr=3e-3, steps=30):
        dtype = torch.float32
        good_rows = np.logical_not(np.isnan(self.sol.a[:,0]))
        good_cols = np.logical_not(np.isnan(self.sol.b[0]))
        u = nn.Parameter(torch.tensor(self.sol.u[:,good_rows],dtype=dtype))
        v = nn.Parameter(torch.tensor(self.sol.v[:,good_cols],dtype=dtype))
        a = nn.Parameter(torch.tensor(self.sol.a[good_rows],dtype=dtype))
        b = nn.Parameter(torch.tensor(self.sol.b[:,good_cols],dtype=dtype))
        o = nn.Parameter(torch.tensor(self.sol.o[:,good_cols],dtype=dtype))
        data = torch.tensor(self.data[good_rows][:,good_cols],dtype=dtype)
        optimizer = torch.optim.Adam([u,v,a,b,o],lr=lr)

        #compute_estimate = lambda u,v,a,b,o : (-2*u.T@v + a + b)**0.5 + o # Had problem with nan values spreading when using this loss space
        compute_part_estimate = lambda u,v,a,b : (-2*u.T@v + a + b)
        
        huberloss = torch.nn.HuberLoss(delta=1)
        for _ in range(steps):

            est = compute_part_estimate(u,v,a,b)
            good_idx = (est ** 0.5 + o).isnan().logical_not()
            normer = 2*np.abs(data - o.detach()) + self.tol**2

            loss = huberloss((est/normer)[good_idx],(((data - o)**2)/normer)[good_idx]) 
            #print("---")
            #print(loss)
            loss += torch.maximum(torch.tensor(0),-est[est.isnan().logical_not()]).mean()
            #print(loss)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        
        self.sol.u[:,good_rows] = u.detach().numpy().astype(dtype=self.sol.u.dtype)
        self.sol.v[:,good_cols] = v.detach().numpy().astype(dtype=self.sol.v.dtype)
        self.sol.a[good_rows] = a.detach().numpy().astype(dtype=self.sol.a.dtype)
        self.sol.b[:,good_cols] = b.detach().numpy().astype(dtype=self.sol.b.dtype)
        self.sol.o[:,good_cols] = o.detach().numpy().astype(dtype=self.sol.o.dtype)

    def upgrade(self, ransac_iter=100): # TODO : is it possible to change what row we do the compaction on?? should be
        u = self.sol.u
        v = self.sol.v
        a = self.sol.a
        b = self.sol.b
        
        u0 = u[:,self.compaction_row:self.compaction_row+1]
        v0 = v[:,self.compaction_col:self.compaction_col+1]
        a = a - 2*u.T@v0 + u0.T@v0
        b = b - 2*u0.T@v + u0.T@v0
        u -= u0
        v -= v0

        diff = (b[:, self.compaction_col] - a[self.compaction_row])
        
        a += diff/2
        b -= diff/2
        c = b[:, self.compaction_col]

        rows_needed = 9
        # B = ((self.data[:, self.compaction_col:self.compaction_col+1] - self.sol.o[:,self.compaction_col])**2 
        #      - (self.data[self.compaction_row, self.compaction_col] - self.sol.o[:,self.compaction_col])**2)
        B = a - c
 
        M = np.stack([u[0,:]**2, 
              2*u[0,:]*u[1,:], 
              2*u[0,:]*u[2,:],
              u[1,:]**2, 
              2*u[1,:]*u[2,:],
              u[2,:]**2, 
              -2*u[0,:],
              -2*u[1,:],
              -2*u[2,:]],axis=1)
        
        most_inliers = -1
        best_sol = None
        for _ in range(ransac_iter):
            temp = np.random.permutation(self.sol.a.shape[0])
            mic_choice = temp[temp != self.compaction_row][:rows_needed]

            sol = np.linalg.solve(M[mic_choice],B[mic_choice])

            H = np.array([[sol[0,0],sol[1,0],sol[2,0]],
                    [sol[1,0],sol[3,0],sol[4,0]],
                    [sol[2,0],sol[4,0],sol[5,0]]
                    ])
            q = sol[6:]
            #u,s,v = np.linalg.svd(np.linalg.inv(H))
            H_inv = np.linalg.inv(H)
            # if not np.all(np.linalg.eigvals(H_inv) > 0):
            #     continue
            # L = np.linalg.cholesky(H_inv).T
            L = np.linalg.cholesky(H_inv + np.eye(3)*(-np.linalg.eigvals(H_inv).min() + 1e-5 if np.linalg.eigvals(H_inv).min() <= 0 else 0)).T
            S = L@v + L@q
            R = np.linalg.solve(L.T,u)

            # L = np.linalg.cholesky(H + np.eye(3)*(-np.linalg.eigvals(H).min() + 1e-5 if np.linalg.eigvals(H).min() <= 0 else 0)).T
            # S = np.linalg.solve(L,v) + q
            # R = L.T@u

            est = np.sqrt(np.sum((np.expand_dims(R.T,1) - np.expand_dims(S.T,0))**2,axis=2)) + self.sol.o
            res = est - self.data
            n_inliers =np.sum(np.abs(res) < self.tol*self.tol_std)
            if n_inliers > most_inliers:
                most_inliers = n_inliers
                best_sol = (R,S)
        print(f"{most_inliers=}, out of {self.data.shape[0]*self.data.shape[1]}")
        if most_inliers == -1:
            raise Exception("Did not find positive definite matrix H in all ransac tries")
        return best_sol

class OffsetSolver(ABC):

    @abstractmethod
    def get_needed_receivers():
        pass

    @abstractmethod
    def get_needed_senders():
        pass

    @abstractmethod
    def solve(data):
        pass
    
class OffsetSolver95(OffsetSolver):

    def get_needed_receivers():
        return 9
    
    def get_needed_senders():
        return 5
    
    def solve(data):
        zsquared = data ** 2
        A = np.concatenate([zsquared[:,1:] - zsquared[:,0:1], -2*data[:,1:], 2*data[:,0:1]],axis=1)
        u = np.linalg.solve(A, np.ones(9))
        sols = np.concatenate([u[-1:]/np.sum(u[:4]),u[4:-1]/u[:4]],axis=0)
        return sols

