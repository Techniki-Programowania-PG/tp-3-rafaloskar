#include <windows.h>
#include <gdiplus.h>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <cmath>
#include <algorithm>

#pragma comment(lib, "gdiplus.lib")

using namespace Gdiplus;

constexpr int numFloors = 5;
constexpr int floorHeight = 100;
constexpr int elevatorOffset = 8;
constexpr int elevatorWidth = 136;
constexpr int elevatorHeight = 90;
constexpr int shaftX = 400;
constexpr int shaftWidth = 150;
constexpr int shaftHeight = numFloors * floorHeight + 70;
constexpr int averagePassengerWeight = 70;
constexpr int maxElevatorLoad = 600;



float elevatorY = (numFloors - 1) * floorHeight + 10;
int currentFloor = 1;
bool moving = false;
// Dodane zmienne globalne:
int doorFrame = 0;
bool doorOpening = false;
bool doorsVisible = false;  // true gdy rysujemy drzwi
int shaftFrame = 0;
bool shaftOpening = false;
bool shaftVisible = false;


struct Passenger {
    int from;
    int to;
    bool inElevator = false;
    bool inElevatorSide = false; // true = prawa
};

struct PassengerExitAnimation {
    int floor;
    bool onRight;
    int frame = 0;
};

struct PassengerEntryAnimation {
    int floor;
    bool fromRight;
    int frame = 0;
};

std::vector<Passenger> passengers;
std::vector<PassengerExitAnimation> animations;
std::vector<PassengerEntryAnimation> entryAnimations;


RECT floorButtons[numFloors][numFloors];
ULONG_PTR gdiplusToken;

int CurrentElevatorLoad() {
    int count = 0;
    for (const auto& p : passengers) {
        if (p.inElevator && p.to >= 1)
            count++;
    }
    return count * averagePassengerWeight;
}


bool IsElevatorEmpty() {
    for (const auto& p : passengers) {
        if (p.inElevator && p.to >= 1)
            return false;
    }
    return true;
}

bool HasPendingRequests() {
    for (const auto& p : passengers) {
        if ((p.inElevator && p.to >= 1) || (!p.inElevator && p.from >= 1 && p.to >= 1))
            return true;
    }
    return false;
}


void DrawFloors(Graphics& g) {
    Pen floorLine(Color(255, 0, 0, 255), 2);   // linie pięter – niebieskie
    Pen shaftPen(Color(255, 0, 0, 0), 3);      // czarny szyb
    Pen eraser(Color(255, 255, 255, 255), 3);  // "gumka" – biała linia

    // --- Rysowanie podłóg pięter ---
    for (int i = 0; i < numFloors; ++i) {
        int y = 30 + (numFloors - i - 1) * floorHeight + elevatorHeight - 20;
        if ((i + 1) % 2 == 1) {
            g.DrawLine(&floorLine, shaftX - 200, y, shaftX, y); // lewa
        }
        else {
            g.DrawLine(&floorLine, shaftX + shaftWidth, y, shaftX + shaftWidth + 200, y); // prawa
        }
    }

    // --- Szyb: góra, dół i stałe ściany ---
    int shaftTopY = 0;
    int shaftBottomY = shaftTopY + shaftHeight + 20;

    g.DrawLine(&shaftPen, shaftX, shaftTopY, shaftX + shaftWidth, shaftTopY);  // góra
    g.DrawLine(&shaftPen, shaftX, shaftBottomY, shaftX + shaftWidth, shaftBottomY);  // dół
    g.DrawLine(&shaftPen, shaftX, shaftTopY, shaftX, shaftBottomY);  // lewa
    g.DrawLine(&shaftPen, shaftX + shaftWidth, shaftTopY, shaftX + shaftWidth, shaftBottomY); // prawa


    // --- Efekt "znikania" tylko po stronie zgodnej z piętrem ---
    // --- Efekt "znikania" i "pojawiania się" jednej ściany szybu ---
    if (shaftVisible) {
        int yTop = static_cast<int>(elevatorY);
        int yBottom = yTop + elevatorHeight;

        bool openRight = (currentFloor % 2 == 0);
        int shaftEdgeX = openRight ? (shaftX + shaftWidth) : shaftX;

        // Animacja otwierania/zamykania
        if (shaftOpening) {
            // OTWIERANIE: biała linia rośnie od dołu do góry
            int yStart = yBottom - shaftFrame;
            g.DrawLine(&eraser, shaftEdgeX, yStart, shaftEdgeX, yBottom);
        }
        else {
            // ZAMYKANIE: biała linia znika od góry do dołu
            int yEnd = yBottom - shaftFrame;
            g.DrawLine(&eraser, shaftEdgeX, yEnd, shaftEdgeX, yBottom);
        }

        // Aktualizacja klatki
        if (shaftOpening && shaftFrame < elevatorHeight) {
            shaftFrame += 2;
        }
        else if (!shaftOpening && shaftFrame > 0) {
            shaftFrame -= 2;
        }
    }
    // --- Prostokąt z obciążeniem windy ---
    int mass = CurrentElevatorLoad();
    int boxWidth = 210;
    int boxHeight = 45;
    int boxX = shaftX + shaftWidth + 40; // po prawej stronie szybu, z odstępem
    int boxY = 30; // wysokość od góry ekranu

    // Tło prostokąta (jasny)
    SolidBrush boxBg(Color(255, 245, 245, 245));
    g.FillRectangle(&boxBg, Rect(boxX, boxY, boxWidth, boxHeight));

    // Obramowanie
    Pen boxPen(Color(255, 80, 80, 80), 2);
    g.DrawRectangle(&boxPen, Rect(boxX, boxY, boxWidth, boxHeight));

    // Tekst
    FontFamily fontFamily(L"Arial");
    Font font(&fontFamily, 16, FontStyleBold);
    SolidBrush textBrush(Color(255, 40, 40, 40));

    wchar_t text[64];
    swprintf(text, 64, L"Obciążenie: %d kg", mass);

    StringFormat format;
    format.SetAlignment(StringAlignmentCenter);     // Wyśrodkowanie poziome
    format.SetLineAlignment(StringAlignmentCenter); // Wyśrodkowanie pionowe

    g.DrawString(
        text, -1, &font,
        RectF((REAL)boxX, (REAL)boxY, (REAL)boxWidth, (REAL)boxHeight),
        &format, &textBrush
    );


}

void DrawElevator(Graphics& g) {

    int x = shaftX + elevatorOffset;
    int y = static_cast<int>(elevatorY);

    Pen pen(Color(255, 255, 0, 0), 2);  // czerwona winda

    // Góra, dół, tylna ściana – zawsze rysowane
    g.DrawLine(&pen, x, y, x + elevatorWidth, y);                         // góra
    g.DrawLine(&pen, x, y + elevatorHeight, x + elevatorWidth, y + elevatorHeight); // dół

    // Drzwi (czyli prawa lub lewa ściana)
    bool openRight = (currentFloor % 2 == 0);  // 2 i 4 → prawa strona otwierana
    int doorX = openRight ? x + elevatorWidth : x;

    int visibleHeight = elevatorHeight - doorFrame;
    if (visibleHeight > 0) {
        int topY = y;
        int bottomY = y + visibleHeight;
        g.DrawLine(&pen, doorX, topY, doorX, bottomY);  // skracana ściana
    }

    // Stała ściana po przeciwnej stronie (nigdy nie znika)
    int otherX = openRight ? x : x + elevatorWidth;
    g.DrawLine(&pen, otherX, y, otherX, y + elevatorHeight);

    // Animacja drzwi
    if (doorsVisible) {
        if (doorOpening && doorFrame < elevatorHeight) {
            doorFrame += 2;
        }
        else if (!doorOpening && doorFrame > 0) {
            doorFrame -= 2;
        }
    }

}

void DrawPassengers(Graphics& g) {
    Pen pen(Color(0, 0, 0));
    SolidBrush brush(Color(0, 0, 0));
    FontFamily fontFamily(L"Arial");
    Font font(&fontFamily, 12);
    SolidBrush textBrush(Color(0, 0, 0));

    auto DrawStickFigure = [&](int x, int y) {
        // Head
        g.FillEllipse(&brush, x, y, 10, 10);
        // Body
        g.DrawLine(&pen, x + 5, y + 10, x + 5, y + 20);
        // Arms
        g.DrawLine(&pen, x, y + 15, x + 10, y + 15);
        // Legs
        g.DrawLine(&pen, x + 5, y + 20, x, y + 30);
        g.DrawLine(&pen, x + 5, y + 20, x + 10, y + 30);
        };

    // --- Czekający pasażerowie przy windzie ---
    for (int floor = 1; floor <= numFloors; ++floor) {
        bool rightSide = (floor % 2 == 0);
        int baseX = rightSide ? (shaftX + shaftWidth + 20) : (shaftX - 30);
        int y = (numFloors - floor) * floorHeight + 60;

        int offset = 0;
        for (const auto& p : passengers) {
            if (!p.inElevator && p.from == floor && p.to >= 1) {
                // Sprawdź czy ten pasażer JEST w trakcie animacji wejścia (czyli jest w entryAnimations)
                bool animEntering = false;
                for (const auto& anim : entryAnimations) {
                    // sprawdzamy po numerze piętra i kierunku (lub można dorzucić docelowe piętro dla unikalności)
                    if (anim.floor == floor && anim.fromRight == (floor % 2 == 0)) {
                        animEntering = true;
                        break;
                    }
                }
                if (animEntering)
                    continue; // NIE rysujemy – jest w trakcie animacji!

                int x = baseX + offset * (rightSide ? 14 : -14);
                DrawStickFigure(x, y);

                std::wstring txt = std::to_wstring(p.to);
                g.DrawString(txt.c_str(), -1, &font,
                    PointF((REAL)x - 4, (REAL)(y - 20)), &textBrush);

                offset++;
            }
        }
    }

    // --- Animacja wchodzenia pasażerów ---
    for (const auto& anim : entryAnimations) {
        int floor = anim.floor;
        bool rightSide = anim.fromRight;
        int xStart = rightSide ? (shaftX + shaftWidth + 10) : (shaftX - 20);
        int yStart = (numFloors - floor) * floorHeight + 60;

        // Docelowe miejsce w windzie – wyśrodkowane z aktualnymi pasażerami
        // Ile osób jest już w windzie (przed tym wchodzącym)?
        int nIn = 0;
        for (const auto& p : passengers)
            if (p.inElevator && p.from == floor)
                nIn++;
        int nAll = nIn + 1; // będzie ten wchodzący

        int figureSpacing = 14;
        int figuresWidth = (nAll > 0) ? (nAll - 1) * figureSpacing : 0;
        int baseElevatorX = shaftX + elevatorOffset + (elevatorWidth - figuresWidth) / 2;
        int targetX = baseElevatorX + nIn * figureSpacing; // pozycja docelowa dla wchodzącego
        int targetY = static_cast<int>(elevatorY) + elevatorHeight - 30;

        float t = anim.frame / 20.0f; // postęp animacji od 0 do 1
        int x = xStart + (int)((targetX - xStart) * t);
        int y = yStart + (int)((targetY - yStart) * t);

        DrawStickFigure(x, y);
    }

    // --- Pasażerowie wewnątrz windy ---
    std::vector<const Passenger*> inside;
    for (const auto& p : passengers) {
        if (p.inElevator && p.to >= 1)
            inside.push_back(&p);
    }
    int n = (int)inside.size();
    int figureSpacing = 14;
    int figuresWidth = (n > 0) ? (n - 1) * figureSpacing : 0;
    int baseX = shaftX + elevatorOffset + (elevatorWidth - figuresWidth) / 2;
    int yIn = static_cast<int>(elevatorY) + elevatorHeight - 30;

    for (int i = 0; i < n; ++i) {
        int x = baseX + i * figureSpacing;
        DrawStickFigure(x, yIn);
        std::wstring txt = std::to_wstring(inside[i]->to);
        g.DrawString(txt.c_str(), -1, &font, PointF((REAL)x - 4, (REAL)(yIn - 20)), &textBrush);
    }

    // --- Animacja wychodzenia pasażerów ---
    for (auto& anim : animations) {
        int y = (numFloors - anim.floor) * floorHeight + 60;
        bool rightSide = (anim.floor % 2 == 0);  // 2 i 4 w prawo, reszta w lewo
        int xStart = rightSide ? shaftX + shaftWidth + 10 : shaftX - 20;
        int x = xStart + anim.frame * (rightSide ? 1 : -1) * 4;
        DrawStickFigure(x, y);
    }

    // Postęp animacji wyjścia
    animations.erase(
        std::remove_if(animations.begin(), animations.end(),
            [](PassengerExitAnimation& a) {
                a.frame++;
                return a.frame > 20;
            }),
        animations.end()
    );
}

void DrawButtons(Graphics& g) {
    FontFamily fontFamily(L"Arial");
    Font font(&fontFamily, 9);
    SolidBrush brush(Color(255, 0, 0, 0)); // 🔴 CZARNY tekst
    Pen border(Color(255, 0, 0, 0));       // Czarna ramka

    for (int i = 0; i < numFloors; ++i) {
        int y = (numFloors - i - 1) * floorHeight + 30;
        int buttonIndex = 0;

        for (int j = 0; j < numFloors; ++j) {
            if (j == i) continue;

            int x = (i % 2 == 0)
                ? (10 + buttonIndex * 35)
                : (850 - (numFloors - 1 - buttonIndex) * 35);

            floorButtons[i][j] = { x, y, x + 30, y + 25 };
            g.DrawRectangle(&border, Rect(x, y, 30, 25));

            std::wstring label = L"  " + std::to_wstring(j + 1);
            g.DrawString(label.c_str(), -1, &font,
                PointF((REAL)(x + 2), (REAL)(y + 4)), &brush);

            buttonIndex++;
        }
    }
}
void CloseDoorsAndShaft(HWND hwnd) {
    doorOpening = false;
    shaftOpening = false;

    while (doorFrame > 0 || shaftFrame > 0) {
        if (doorFrame > 0) doorFrame--;
        if (shaftFrame > 0) shaftFrame--;
        InvalidateRect(hwnd, NULL, FALSE);
        UpdateWindow(hwnd);
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    doorsVisible = false;
    shaftVisible = false;
}
void MoveElevator(HWND hwnd) {
    if (moving) return;
    moving = true;

    while (true) {
        // Sprawdź, czy są zadania: ktoś w windzie LUB ktoś czeka na piętrze i da się go zabrać
        bool hasTasks = false;
        int load = CurrentElevatorLoad();
        for (const auto& p : passengers) {
            if (p.inElevator && p.to >= 1) {
                hasTasks = true;
                break;
            }
        }
        // Czy można zabrać kogokolwiek z poczekalni?
        for (const auto& p : passengers) {
            if (!p.inElevator && p.from >= 1 && p.to >= 1 && load + averagePassengerWeight <= maxElevatorLoad) {
                hasTasks = true;
                break;
            }
        }
        if (!hasTasks) {
            // Gdy NIE ma zgłoszeń: sprawdź, czy jesteśmy pustą windą na nieparterze
            if (IsElevatorEmpty() && currentFloor != 1) {
                bool interrupted = false;
                for (int i = 0; i < 50; ++i) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    InvalidateRect(hwnd, NULL, FALSE);
                    UpdateWindow(hwnd);
                    // Jeśli pojawi się zadanie lub ktoś wsiadł, przerwij czekanie
                    if (!IsElevatorEmpty() || HasPendingRequests()) {
                        interrupted = true;
                        break;
                    }
                }
                if (interrupted)
                    continue;
                CloseDoorsAndShaft(hwnd);
                float targetY = (numFloors - 1) * floorHeight + 10;
                while (std::abs(elevatorY - targetY) > 1.0f) {
                    elevatorY += (elevatorY < targetY) ? 1.0f : -1.0f;
                    InvalidateRect(hwnd, NULL, FALSE);
                    UpdateWindow(hwnd);
                    std::this_thread::sleep_for(std::chrono::milliseconds(5));
                }
                elevatorY = targetY;
                currentFloor = 1;
            }
            break;
        }

        // Znajdź najbliższy cel:
        int targetFloor = -1;
        int minDist = numFloors + 1;
        // Najpierw do pasażera w windzie
        for (const auto& p : passengers) {
            if (p.inElevator && p.to >= 1) {
                int dist = std::abs(currentFloor - p.to);
                if (dist < minDist) {
                    minDist = dist;
                    targetFloor = p.to;
                }
            }
        }
        // Potem do oczekujących, których możemy zabrać (tylko tylu, ilu zmieści się wagowo!)
        int tmpLoad = load;
        for (const auto& p : passengers) {
            if (!p.inElevator && p.from >= 1 && p.to >= 1 && tmpLoad + averagePassengerWeight <= maxElevatorLoad) {
                int dist = std::abs(currentFloor - p.from);
                if (dist < minDist) {
                    minDist = dist;
                    targetFloor = p.from;
                }
                tmpLoad += averagePassengerWeight; // symulujemy kolejnych pasażerów do zabrania
            }
        }
        if (targetFloor == -1) break;

        // Ruch do piętra
        float targetY = (numFloors - targetFloor) * floorHeight + 10;
        while (std::abs(elevatorY - targetY) > 1.0f) {
            elevatorY += (elevatorY < targetY) ? 1.0f : -1.0f;
            InvalidateRect(hwnd, NULL, FALSE);
            UpdateWindow(hwnd);
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        elevatorY = targetY;
        currentFloor = targetFloor;

        // OTWIERANIE drzwi i tunelu
        doorOpening = true;
        doorsVisible = true;
        shaftOpening = true;
        shaftVisible = true;
        doorFrame = 0;
        shaftFrame = 0;
        while (doorFrame < elevatorHeight && shaftFrame < elevatorHeight) {
            doorFrame++;
            shaftFrame++;
            InvalidateRect(hwnd, NULL, FALSE);
            UpdateWindow(hwnd);
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }

        // Pasażerowie wysiadają
        for (auto& p : passengers) {
            if (p.inElevator && p.to == currentFloor) {
                p.inElevator = false;
                p.to = -1;
                animations.push_back({ currentFloor, p.inElevatorSide, 0 });
            }
        }

        int nLoad = CurrentElevatorLoad();
        for (auto& p : passengers) {
            if (!p.inElevator && p.from == currentFloor && p.to >= 1 && nLoad + averagePassengerWeight <= maxElevatorLoad) {
                // Dodaj animację wejścia
                entryAnimations.push_back({ currentFloor, (currentFloor % 2 == 0), 0 });
                // Animuj wejście
                for (int i = 0; i < 20; ++i) {
                    entryAnimations.back().frame = i;
                    InvalidateRect(hwnd, NULL, FALSE);
                    UpdateWindow(hwnd);
                    std::this_thread::sleep_for(std::chrono::milliseconds(18));
                }
                // Po animacji dopiero "wchodzi" do windy
                p.inElevator = true;
                p.inElevatorSide = (currentFloor % 2 == 0);
                nLoad += averagePassengerWeight;
                entryAnimations.pop_back();
            }
        }



        CloseDoorsAndShaft(hwnd);

        // Usuwanie pasażerów, którzy już wysiedli
        passengers.erase(
            std::remove_if(passengers.begin(), passengers.end(), [](const Passenger& p) {
                return !p.inElevator && (p.to == -1 || p.from == -1);
                }),
            passengers.end()
        );
    }

    moving = false;
    InvalidateRect(hwnd, NULL, FALSE);
}
void AddFloorCall(int from, int to, HWND hwnd) {
    passengers.push_back({ from, to });
    if (!moving) {
        std::thread(MoveElevator, hwnd).detach();
    }
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_LBUTTONDOWN: {
        int x = LOWORD(lParam);
        int y = HIWORD(lParam);

        for (int i = 0; i < numFloors; ++i) {
            for (int j = 0; j < numFloors; ++j) {
                if (j == i) continue;
                RECT r = floorButtons[i][j];
                if (x >= r.left && x <= r.right && y >= r.top && y <= r.bottom) {
                    AddFloorCall(i + 1, j + 1, hwnd);
                }
            }
        }
        break;
    }
    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);

        RECT rect;
        GetClientRect(hwnd, &rect);
        int width = rect.right - rect.left;
        int height = rect.bottom - rect.top;

        // Tworzenie bufora w pamięci
        HDC memDC = CreateCompatibleDC(hdc);
        HBITMAP memBitmap = CreateCompatibleBitmap(hdc, width, height);
        HBITMAP oldBitmap = (HBITMAP)SelectObject(memDC, memBitmap);

        Graphics g(memDC);
        g.Clear(Color(255, 255, 255)); // Białe tło

        DrawFloors(g);
        DrawButtons(g);
        DrawElevator(g);
        DrawPassengers(g);

        // Przekopiowanie z bufora na ekran
        BitBlt(hdc, 0, 0, width, height, memDC, 0, 0, SRCCOPY);

        // Czyszczenie
        SelectObject(memDC, oldBitmap);
        DeleteObject(memBitmap);
        DeleteDC(memDC);

        EndPaint(hwnd, &ps);
        break;
    }
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
    return 0;
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow) {
    GdiplusStartupInput gdiplusStartupInput;
    GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

    const wchar_t CLASS_NAME[] = L"ElevatorSim";
    WNDCLASS wc = {};
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;
    RegisterClass(&wc);

    HWND hwnd = CreateWindowEx(0, CLASS_NAME, L"Symulator windy", WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 1000, 700, NULL, NULL, hInstance, NULL);
    ShowWindow(hwnd, nCmdShow);

    MSG msg = {};
    while (GetMessage(&msg, NULL, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    GdiplusShutdown(gdiplusToken);
    return 0;
}
