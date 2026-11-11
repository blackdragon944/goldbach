"""
Треугольник Серпинского с бесконечным приближением.

Запуск:   python sierpinski.py
Нужно:    pip install matplotlib numpy

Управление в окне:
  колесо мыши        — приближение / отдаление к курсору
  левая кнопка+тянуть — перемещение
  стрелки  ↑ / ↓      — вырастить / уменьшить на один шаг (или ползунок внизу)
  клавиша  r          — вернуть исходный вид
  клавиши  + / -      — приближение / отдаление к центру
  клавиша  s          — сохранить картинку в PNG

Все параметры — в блоке НАСТРОЙКИ ниже.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

# =====================================================================
#                             НАСТРОЙКИ
# =====================================================================

# --- Самый маленький треугольник ------------------------------------
# Меньше него ничего не бывает. Фигура из него только РАСТЁТ.
V1 = (0.0, 3.0)             # вершина 1  (из неё идёт рост, она не двигается)
V2 = (0.0, 4.0)             # вершина 2
V3 = (1.0, 4.0)             # вершина 3

K = 4                       # шаг роста при запуске: 3**K маленьких
                            # треугольников, размер больше в 2**K раз.
                            # В окне меняется ползунком или стрелками ↑ ↓
K_MAX = 20

# --- Поворот и зеркало (применяются ко всей фигуре) -----------------
ROTATE_DEG = 0.0            # угол поворота, против часовой стрелки
PIVOT = 1                   # вокруг какой вершины крутить: 1, 2 или 3
MIRROR = None               # None, "x" (верх-низ), "y" (лево-право),
                            # "y=x" или "y=-x"  (ось проходит через PIVOT)

# --- Вид фрактала -----------------------------------------------------
INVERT = False              # True: закрасить дырки вместо треугольников
SHOW_OUTLINE = True         # пунктир вокруг большого треугольника

# --- Сколько треугольников и где -----------------------------------
# Каждая строка — одна копия: (сдвиг по x, сдвиг по y, заливка, контур).
# Сдвиг считается от треугольника V1, V2, V3. Все копии растут одинаково.
COPIES = [
    (0, 0, "#7F77DD", "#534AB7"),   # 1-й: (0, 3), (0, 4), (1, 2)
    (0, 1, "#5DCAA5", "#0F6E56"),   # 2-й: на 1 выше: (0, 4), (0, 5), (1, 3)
    # (0, 2, "#F0997B", "#993C1D"), # 3-й: ещё на 1 выше — раскомментируйте
]
MIN_PIXELS = 1.5            # до какого размера (в пикселях) дробить
MAX_TRIANGLES = 200_000     # предел на один кадр, чтобы не тормозило

# --- Ваши графики -----------------------------------------------------
# Любые функции y = f(x). Пишите через numpy: np.sin(x), x**2, ...
MY_GRAPHS = [
    (lambda x: x + 3,  "y = x + 3",  "#E24B4A"),
    (lambda x: -x + 3, "y = -x + 3", "#1D9E75"),
    # (lambda x: np.sin(x) + 5, "y = sin x + 5", "#378ADD"),
]
# Вертикальные прямые x = c
MY_VERTICAL = [
    # (0, "x = 0", "#888780"),
]
# Отрезать всё, что правее/ниже какой-то линии.
# Функция возвращает True для точек, которые ОСТАВЛЯЕМ.
# Пример: оставить только левее y = x + 3  ->  lambda x, y: y >= x + 3
KEEP = None
# KEEP = lambda x, y: y >= x + 3 - 1e-12

# --- Начальный вид ------------------------------------------------
VIEW = None                 # None = подогнать под треугольник,
                            # или (xmin, xmax, ymin, ymax)
GRID = True

# =====================================================================
#                      дальше можно не трогать
# =====================================================================


def transform(points):
    """Поворот и зеркало вокруг выбранной вершины."""
    pts = np.array(points, dtype=float)
    O = pts[PIVOT - 1].copy()
    p = pts - O
    a = math.radians(ROTATE_DEG)
    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    p = p @ R.T
    if MIRROR == "x":
        p[:, 1] *= -1
    elif MIRROR == "y":
        p[:, 0] *= -1
    elif MIRROR == "y=x":
        p = p[:, ::-1].copy()
    elif MIRROR == "y=-x":
        p = -p[:, ::-1]
    return p + O


SMALL = transform([V1, V2, V3])          # самый маленький треугольник
ANCHOR = SMALL[0].copy()                 # точка роста = вершина 1


def big_triangle(k, off=(0.0, 0.0)):
    """Маленький треугольник, увеличенный в 2**k раз от вершины 1, + сдвиг."""
    return ANCHOR + (SMALL - ANCHOR) * (2 ** k) + np.array(off, dtype=float)


state = {"k": K}


def tri_bbox(t):
    xs, ys = t[:, 0], t[:, 1]
    return xs.min(), xs.max(), ys.min(), ys.max()


def keep_tri(t):
    if KEEP is None:
        return True
    c = t.mean(axis=0)
    return bool(KEEP(c[0], c[1]))


def collect(view, px_size, off=(0.0, 0.0)):
    """Обходит фрактал и берёт только то, что видно на экране."""
    x0, x1, y0, y1 = view
    tris, holes = [], []
    stack = [(big_triangle(state["k"], off), state["k"])]
    limit = MIN_PIXELS * px_size
    while stack:
        t, d = stack.pop()
        bx0, bx1, by0, by1 = tri_bbox(t)
        if bx1 < x0 or bx0 > x1 or by1 < y0 or by0 > y1:
            continue                       # вне экрана — пропускаем
        size = max(bx1 - bx0, by1 - by0)
        if d == 0 or size <= limit or len(tris) + len(stack) > MAX_TRIANGLES:
            if keep_tri(t):                # d == 0: самый маленький, дальше
                tris.append(t)             # не делим
            continue
        a, b, c = t
        ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
        hole = np.array([ab, bc, ca])
        if INVERT and keep_tri(hole):
            holes.append(hole)
        stack.append((np.array([a, ab, ca]), d - 1))
        stack.append((np.array([ab, b, bc]), d - 1))
        stack.append((np.array([ca, bc, c]), d - 1))
    return tris, holes


from matplotlib.widgets import Slider

fig, ax = plt.subplots(figsize=(8, 8.6))
fig.subplots_adjust(bottom=0.1, top=0.97, left=0.08, right=0.97)
fig.canvas.manager.set_window_title("Треугольник Серпинского")
ax.set_aspect("equal")
if GRID:
    ax.grid(True, color="#dddddd", linewidth=0.6)
ax.axhline(0, color="#999999", linewidth=0.8)
ax.axvline(0, color="#999999", linewidth=0.8)

layers = []      # для каждой копии: (сдвиг, заливка, контур, коллекция, пунктир)
for dx, dy, fc, ec in COPIES:
    coll = PolyCollection([], linewidths=0.3, alpha=0.85)
    ax.add_collection(coll)
    (ol,) = ax.plot([], [], "--", color=ec, linewidth=1, visible=SHOW_OUTLINE)
    layers.append(((dx, dy), fc, ec, coll, ol))
    # вершины самого маленького треугольника этой копии (не двигаются)
    for i, v in enumerate(SMALL + np.array([dx, dy], dtype=float)):
        ax.plot(*v, "o", color="#D85A30" if i == 0 else ec, ms=5, zorder=5)
        ax.annotate(f"({v[0]:g}, {v[1]:g})", v, textcoords="offset points",
                    xytext=(6, 6), fontsize=9, color=ec, zorder=5)

graph_lines = []
for f, label, color in MY_GRAPHS:
    (ln,) = ax.plot([], [], color=color, linewidth=1.5, label=label, zorder=4)
    graph_lines.append((f, ln))
for c, label, color in MY_VERTICAL:
    ax.axvline(c, color=color, linewidth=1.5, label=label, zorder=4)
if MY_GRAPHS or MY_VERTICAL:
    ax.legend(loc="upper right")

info = ax.text(0.01, 0.01, "", transform=ax.transAxes, fontsize=9,
               va="bottom", ha="left",
               bbox=dict(boxstyle="round", fc="white", ec="#cccccc", alpha=0.85))

def home_view():
    if VIEW is not None:
        return VIEW
    allp = np.vstack([big_triangle(state["k"], (dx, dy))
                      for dx, dy, *_ in COPIES])
    bx0, bx1, by0, by1 = tri_bbox(allp)
    w, h = bx1 - bx0, by1 - by0
    s = max(w, h) * 1.16
    cx, cy = (bx0 + bx1) / 2, (by0 + by1) / 2
    return (cx - s / 2, cx + s / 2, cy - s / 2, cy + s / 2)


HOME = home_view()
HOME_W = HOME[1] - HOME[0]

busy = False


def redraw(_=None):
    global busy
    if busy:
        return
    busy = True
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    w_px = max(ax.get_window_extent().width, 1)
    px = (x1 - x0) / w_px
    k = state["k"]
    vert_lines = []
    for n, (off, fc, ec, coll, ol) in enumerate(layers, 1):
        tris, holes = collect((x0, x1, y0, y1), px, off)
        if INVERT:
            coll.set_verts(holes)
            coll.set_facecolor(fc)
            coll.set_edgecolor("none")
        else:
            coll.set_verts(tris)
            coll.set_facecolor(fc)
            coll.set_edgecolor(ec if len(tris) < 3000 else "none")
        B = big_triangle(k, off)
        o = np.vstack([B, B[:1]])
        ol.set_data(o[:, 0], o[:, 1])
        vert_lines.append(f"{n}: " + ";  ".join(f"({p[0]:g}, {p[1]:g})"
                                                for p in B))
    xs = np.linspace(x0, x1, 2000)
    for f, ln in graph_lines:
        with np.errstate(all="ignore"):
            ys = f(xs) * np.ones_like(xs)
        ln.set_data(xs, ys)
    zoom = HOME_W / (x1 - x0)
    txt = (f"шаг {k}: в каждом {3**k:_} маленьких, размер ×{2**k:_}\n"
           + "\n".join(vert_lines) +
           f"\nприближение ×{zoom:_.0f}")
    info.set_text(txt.replace("_", " "))
    fig.canvas.draw_idle()
    busy = False


def set_view(x0, x1, y0, y1):
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    redraw()


def zoom_at(x, y, factor):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if (x1 - x0) * factor < 1e-12 * max(1.0, abs(x), abs(y)):
        print("Дальше не приблизить: кончилась точность чисел.")
        return
    set_view(x - (x - x0) * factor, x + (x1 - x) * factor,
             y - (y - y0) * factor, y + (y1 - y) * factor)


def on_scroll(e):
    if e.inaxes != ax:
        return
    zoom_at(e.xdata, e.ydata, 0.7 if e.button == "up" else 1 / 0.7)


drag = {}


def on_press(e):
    tb = fig.canvas.toolbar
    if e.inaxes == ax and e.button == 1 and (tb is None or tb.mode == ""):
        drag["p"] = (e.x, e.y, ax.get_xlim(), ax.get_ylim())


def on_move(e):
    if "p" not in drag:
        return
    sx, sy, (x0, x1), (y0, y1) = drag["p"]
    bb = ax.get_window_extent()
    dx = (e.x - sx) * (x1 - x0) / bb.width
    dy = (e.y - sy) * (y1 - y0) / bb.height
    ax.set_xlim(x0 - dx, x1 - dx)
    ax.set_ylim(y0 - dy, y1 - dy)
    fig.canvas.draw_idle()


def on_release(e):
    if drag.pop("p", None) is not None:
        redraw()


def set_k(k, refit=True):
    """Новый шаг роста. refit=True — подогнать вид под новый размер."""
    global HOME, HOME_W
    k = int(max(0, min(K_MAX, k)))
    state["k"] = k
    HOME = home_view()
    HOME_W = HOME[1] - HOME[0]
    if int(k_slider.val) != k:
        k_slider.eventson = False
        k_slider.set_val(k)
        k_slider.eventson = True
    if refit:
        set_view(*HOME)
    else:
        redraw()


sax = fig.add_axes([0.18, 0.025, 0.62, 0.03])
k_slider = Slider(sax, "шаг роста", 0, K_MAX, valinit=K, valstep=1,
                  color="#7F77DD")
k_slider.on_changed(lambda v: set_k(v) if int(v) != state["k"] else None)


def on_key(e):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    if e.key == "r":
        set_view(*HOME)
    elif e.key == "up":
        set_k(state["k"] + 1)
    elif e.key == "down":
        set_k(state["k"] - 1)
    elif e.key in ("+", "="):
        zoom_at(cx, cy, 0.5)
    elif e.key == "-":
        zoom_at(cx, cy, 2)
    elif e.key == "s":
        name = "sierpinski_view.png"
        fig.savefig(name, dpi=200)
        print("Сохранено:", name)


fig.canvas.mpl_connect("scroll_event", on_scroll)
fig.canvas.mpl_connect("button_press_event", on_press)
fig.canvas.mpl_connect("motion_notify_event", on_move)
fig.canvas.mpl_connect("button_release_event", on_release)
fig.canvas.mpl_connect("key_press_event", on_key)
fig.canvas.mpl_connect("resize_event", redraw)
ax.ticklabel_format(useOffset=False)
# лупа на панели инструментов тоже работает: после неё перерисовываем
fig.canvas.mpl_connect("button_release_event",
                       lambda e: redraw() if fig.canvas.toolbar
                       and fig.canvas.toolbar.mode != "" else None)

set_view(*HOME)
plt.show()
