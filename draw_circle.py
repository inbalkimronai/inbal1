import turtle


def main():
    try:
        diameter = float(input("הזן קוטר לעיגול: "))
        if diameter <= 0:
            print("הקוטר חייב להיות מספר חיובי.")
            return
    except ValueError:
        print("נא הזן מספר תקין.")
        return

    radius = diameter / 2

    screen = turtle.Screen()
    screen.title("עיגול לפי קוטר")
    screen.bgcolor("white")

    pen = turtle.Turtle()
    pen.speed(0)
    pen.color("blue")
    pen.pensize(2)

    pen.penup()
    pen.goto(0, -radius)
    pen.pendown()
    pen.circle(radius)

    pen.hideturtle()
    screen.mainloop()


if __name__ == "__main__":
    main()
