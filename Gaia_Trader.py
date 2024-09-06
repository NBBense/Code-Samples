import numpy as np
import pandas as pd
import yfinance as yf
from scipy.optimize import minimize
from ib_insync import *

# Define the stock symbols and the time period
stocks = ['AAPL', 'MSFT', 'GOOGL', 'AMZN', 'TSLA']
start_date = '2022-01-01'
end_date = '2023-01-01'

# Download stock data
data = yf.download(stocks, start=start_date, end=end_date)['Adj Close']
returns = data.pct_change().dropna()

# Define the objective function (negative Sharpe ratio)
def objective_function(weights):
    portfolio_return = np.sum(returns.mean() * weights) * 252
    portfolio_volatility = np.sqrt(np.dot(weights.T, np.dot(returns.cov() * 252, weights)))
    sharpe_ratio = portfolio_return / portfolio_volatility
    return -sharpe_ratio

# Define genetic algorithm parameters
population_size = 100
num_generations = 50
mutation_rate = 0.02

# Initialize population
def initialize_population(size, num_stocks):
    return np.random.dirichlet(np.ones(num_stocks), size=size)

# Crossover operation
def crossover(parent1, parent2):
    point = np.random.randint(1, len(parent1) - 1)
    child1 = np.concatenate([parent1[:point], parent2[point:]])
    child2 = np.concatenate([parent2[:point], parent1[point:]])
    return child1, child2

# Mutation operation
def mutate(individual, rate):
    mutation = np.random.rand(len(individual)) < rate
    individual[mutation] = np.random.rand(np.sum(mutation))
    individual /= np.sum(individual)
    return individual

# Main genetic algorithm loop
def genetic_algorithm():
    num_stocks = len(stocks)
    population = initialize_population(population_size, num_stocks)
    for generation in range(num_generations):
        # Evaluate fitness
        fitness = np.array([objective_function(ind) for ind in population])
        sorted_indices = np.argsort(fitness)
        population = population[sorted_indices]
        
        # Select the best individuals
        elite = population[:10]
        
        # Generate new population
        new_population = []
        while len(new_population) < population_size:
            parents = np.random.choice(elite, size=2, replace=False)
            child1, child2 = crossover(parents[0], parents[1])
            new_population.append(mutate(child1, mutation_rate))
            new_population.append(mutate(child2, mutation_rate))
        
        population = np.array(new_population)
        print(f'Generation {generation + 1} - Best Fitness: {-fitness[sorted_indices][0]}')

    best_individual = population[np.argmax(fitness)]
    return best_individual

# Run the genetic algorithm
best_portfolio_weights = genetic_algorithm()
print('Best Portfolio Weights:', best_portfolio_weights)

# Portfolio performance metrics
portfolio_return = np.sum(returns.mean() * best_portfolio_weights) * 252
portfolio_volatility = np.sqrt(np.dot(best_portfolio_weights.T, np.dot(returns.cov() * 252, best_portfolio_weights)))
sharpe_ratio = portfolio_return / portfolio_volatility
print(f'Expected Annual Return: {portfolio_return:.2f}')
print(f'Expected Annual Volatility: {portfolio_volatility:.2f}')
print(f'Sharpe Ratio: {sharpe_ratio:.2f}')

# Connect to IB
ib = IB()
ib.connect('127.0.0.1', 7497, 1)  # Replace with IB Gateway/TWS connection details

# Define the stock contracts
def get_stock_contracts(stocks):
    contracts = []
    for stock in stocks:
        contract = Stock(stock, 'SMART', 'USD')
        contracts.append(contract)
    return contracts

contracts = get_stock_contracts(stocks)

# Place trades based on the optimized portfolio
def place_trades(weights, contracts):
    for weight, contract in zip(weights, contracts):
        if weight > 0:
            order = MarketOrder('BUY', 100,  # Example: buy 100 shares
                                orderType='MKT')
            trade = ib.placeOrder(contract, order)
            print(f'Placed order: {weight*100}% of {contract.symbol}')
        else:
            print(f'Skipping {contract.symbol} due to zero weight.')

place_trades(best_portfolio_weights, contracts)

# Disconnect from IB
ib.disconnect()
