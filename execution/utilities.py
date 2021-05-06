def banner(s, width=69):
    stars = '*' * width
    pad = (width + len(s)) // 2
    print(f'{stars}\n{""}\n{s:>{pad}}\n{""}\n{stars}')


def writeVersionHeader(file):
    version = "pyFSI Output File - Martin Saravia - 2021"
    head = "#" + (67 * "-") + "#\n"
    pad = ((67 - len(version)) // 2) * " "
    message = "#" + pad + version + pad + "#"
    file.write(head)
    file.write(message + "\n")
    file.write(head)

