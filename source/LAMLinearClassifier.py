import torch
import numpy as np
import torch.nn.functional as F


class My_model(torch.nn.Module):
    def __init__(self, input_dimension):
        super().__init__()
        self.linear = torch.nn.Linear(input_dimension, 1)

    def forward(self, input_dimension):
        hidden = F.sigmoid(self.linear(input_dimension))
        hidden_1 = 1 - hidden
        hidden_2 = hidden - 0
        output = torch.cat([hidden_1, hidden_2], dim=1)
        return output


def criterion(output_, y):  # output: 10000 * 2 soft logits
    # print(output[0])
    num = output_.size(0)
    output = F.gumbel_softmax(output_, tau=1, hard=True) # output: 10000 * 2 one-hot [1,0] [0,1]
    # print(output[0])
    all_res = []
    for i in range(num):
        y1 = torch.Tensor([y[i] ** 2, -y[i] ** 2])
        y2 = torch.Tensor([-y[i] ** 2, y[i] ** 2])
        output_item = output[i]
        if y[i].numpy().tolist()[0] > 0:
            res = torch.sum(output_item * y1).unsqueeze(0)
        else:
            res = torch.sum(output_item * y2).unsqueeze(0)
        all_res.append(res)
    all_res = torch.cat(all_res, dim=0)
    loss = torch.mean(all_res)

    gumble = sum(np.argmax(output.detach().numpy(), axis=1))
    return loss, gumble


def configure_optimizer(model):
    return torch.optim.Adam(model.parameters(), lr=1e-2)


def full_gd(model, optimizer, X_train, y_train, n_epochs=2000):
    train_losses = np.zeros(n_epochs)

    for it in range(n_epochs):
        outputs = model(X_train)
        loss, gumble = criterion(outputs, y_train)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        train_losses[it] = loss.item()
    return train_losses


# torch.manual_seed(2022)
x = []
y = []
with open("gamma/youtube100_0.9.lamfinetune", 'r') as f:
    line = f.readline()
    while line:
        line = line.split(' ')
        x.append([])
        x[len(x)-1].append(float(line[0]))
        x[len(x)-1].append(float(line[1]))
        y.append(float(line[2]))
        line =f.readline()

x_np = np.array(x)
y_np = np.array(y)

_, input_dimension = x_np.shape

model = My_model(input_dimension)

X_train = torch.from_numpy(x_np.astype(np.float32))
y_train = torch.from_numpy(y_np.astype(np.float32)).reshape(-1, 1)

optimizer = configure_optimizer(model)
train_losses = full_gd(model, optimizer, X_train, y_train)

# plt.plot(train_losses, label='train loss')
# plt.legend()
# plt.show()

print("parameters:")
print(model.linear.weight)
print(model.linear.bias)
