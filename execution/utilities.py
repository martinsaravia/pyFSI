def banner(s, width=69):
    stars = '*' * width
    pad = (width + len(s)) // 2
    print(f'{stars}\n{""}\n{s:>{pad}}\n{""}\n{stars}')
