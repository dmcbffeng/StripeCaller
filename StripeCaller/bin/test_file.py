#!/usr/bin/env python

import argparse


def MyParser():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        '--s'
    )
    args, _ = parse.parse_known_args()
    return args


if __name__ == '__main__':
    print('test')
    p = MyParser()
    print(p.s)

