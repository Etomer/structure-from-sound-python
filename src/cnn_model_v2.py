import torch
import torch.nn as nn
import lightning as L
from src.model_pieces import Block


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

class cnn_model_v2(L.LightningModule):

    def __init__(self,
         n_shift_bins: int = 500,
         dropout: float = 0.1,
         learning_rate: float = 3e-4,
         scale_cnn_width: int = 1,
         n_frequency_components_input: int = 2500,
         loss_fn = nn.CrossEntropyLoss(),
         n_blocks : int = 2,
         block_width : int = 1000,
         ):
        super().__init__()
        self.save_hyperparameters(ignore=['loss_fn'])
        self.learning_rate = learning_rate
        self.loss_fn  = loss_fn
        
        self.width_at_scale_1 = 576 # Could compute this value with formula depending on self.cnn bellow. Easier just to test 

        self.flatten_and_project_to_block_dim = nn.Sequential(
            nn.Flatten(),
            nn.Dropout(dropout),
            nn.Linear(scale_cnn_width*self.width_at_scale_1,block_width),
            nn.GELU(),
        )
        self.blocks = nn.Sequential(*[Block(block_width,dropout) for i in range(n_blocks)])
        self.classifier = nn.Linear(block_width,n_shift_bins)
        
        self.apply(self._init_weights)
        
        self.cnn = nn.Sequential(
            nn.Conv1d(4,48*scale_cnn_width, 50,stride=5),
            nn.GELU(),
            nn.Conv1d(48*scale_cnn_width,48*scale_cnn_width, 50,stride=5),
            nn.GELU(),
            nn.Conv1d(48*scale_cnn_width,48*scale_cnn_width, 30,stride=5),
            nn.GELU(),
            nn.Flatten(),
        )
        

    def forward(self, x):
        x /= x.std(dim=(1,2),keepdim=True) + 1e-5
        x = self.cnn(x)
        x = self.flatten_and_project_to_block_dim(x)
        x = self.blocks(x)
        x = self.classifier(x)
        return x

    def training_step(self, batch):
        x, y = batch
        logits = self(x)
        loss = self.loss_fn(logits, y)
        self.log('train/loss', loss, on_epoch=True)
        return loss
    
    def validation_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x)
        loss = self.loss_fn(logits, y)
        self.log("val/loss_epoch", loss, on_step=False, on_epoch=True)
        pass

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.parameters(), lr=self.learning_rate)
        return optimizer


    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.0002)
            if module.bias is not None:
                torch.nn.init.zeros_(module.bias)
        elif isinstance(module, nn.Embedding):
            torch.nn.init.normal_(module.weight, mean=0.0, std=0.0002)
        
    