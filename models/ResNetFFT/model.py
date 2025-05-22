import torch
import torch.nn as nn

class ResNetFFT(nn.Module):

    def __init__(self, input_size, block_width, n_blocks, output_size, dropout = 0.1, scale_cnn = 1):
        super().__init__()
        self.cnn = nn.Sequential(
            nn.Conv1d(4,10*scale_cnn,50,stride=5),
            nn.GELU(),
            nn.Conv1d(10*scale_cnn,30*scale_cnn,50,stride=5),
            nn.GELU(),
            nn.Conv1d(30*scale_cnn,50*scale_cnn,30,stride=5),
            nn.GELU(),
            nn.Flatten(),
        )
        
        self.proj = nn.Linear(800*scale_cnn,block_width)
        self.thinker = nn.Sequential(*[Block(block_width,dropout) for i in range(n_blocks)])
        self.classifier = nn.Linear(block_width,output_size)
        self.apply(self._init_weights)

    def forward(self, x):
        temp = torch.fft.fft(x)
        temp = temp[:,:,:3000]
        x = torch.concatenate([temp.real,temp.imag],dim=1)
        x /= x.std(dim=(1,2),keepdim=True) + 1e-5
        x = self.cnn(x)
        x = self.proj(x)
        x = self.thinker(x)
        x = self.classifier(x)
        return x
        
    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.0002)
            if module.bias is not None:
                torch.nn.init.zeros_(module.bias)
        elif isinstance(module, nn.Embedding):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.0002)

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
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.0002)
            if module.bias is not None:
                torch.nn.init.zeros_(module.bias)
        elif isinstance(module, nn.Embedding):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.0002)
            
    def forward(self, x):
        return x + self.l2(self.act(self.l(self.ln(self.dropout(x)))))