import torch
import torch.nn as nn

class Block(nn.Module):
    def __init__(self,size, dropout):
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        self.l = nn.Linear(size,2*size)
        self.l2 = nn.Linear(2*size,size)
        self.act = nn.GELU()
        self.ln = nn.LayerNorm(size)
        self.apply(self._init_weights)
        
    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.002)
            if module.bias is not None:
                torch.nn.init.zeros_(module.bias)
        elif isinstance(module, nn.Embedding):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.002)
            
    def forward(self, x):
        return x + self.l2(self.act(self.l(self.ln(self.dropout(x)))))