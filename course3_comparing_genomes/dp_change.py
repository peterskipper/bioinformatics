from fire import Fire
import numpy as np


def dp_change(money, coin_arr):
    MIN_NUM_COINS = {0: 0}
    coin_arr = sorted(coin_arr, reverse=True)
    for amt in range(1, money + 1):
        min_num = np.inf
        for coin in coin_arr:
            if amt >= coin:
                if MIN_NUM_COINS[amt - coin] + 1 < min_num:
                    min_num = MIN_NUM_COINS[amt - coin] + 1
        MIN_NUM_COINS[amt] = min_num
    return MIN_NUM_COINS[money]


def main(money, coin_list):
    money = int(money)
    coin_arr = [int(coin) for coin in coin_list]
    return dp_change(money, coin_arr)


if __name__ == '__main__':
    Fire(main)
